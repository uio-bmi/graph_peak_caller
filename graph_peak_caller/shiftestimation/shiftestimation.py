# Time-stamp: <2015-03-10 15:38:17 Tao Liu>

"""Module Description
Copyright (c) 2008,2009 Yong Zhang, Tao Liu <taoliu@jimmy.harvard.edu>
Copyright (c) 2010,2011 Tao Liu <taoliu@jimmy.harvard.edu>
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).
@status:  experimental
@version: $Revision$
@author:  Yong Zhang, Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""
import numpy as np
import logging
from collections import defaultdict

def median(nums):
    """Calculate Median.
    Parameters:
    nums:  list of numbers
    Return Value:
    median value
    """
    p = sorted(nums)
    l = len(p)
    if l % 2 == 0:
        return (p[l//2]+p[l//2-1])//2
    else:
        return p[l//2]


class NotEnoughPairsException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Treatment:
    def __init__(self, dicts):
        assert isinstance(dicts, dict)
        assert "+" in dicts and "-" in dicts and len(dicts) == 2

        self._chroms = list(set(list(dicts["+"].keys()) + list(dicts["-"].keys())))
        logging.info("Treatment init")
        logging.info("Chroms: %s" % self._chroms)
        self._pos_dict = {key: np.array(dicts["+"][key], dtype="int")
                          for key in self._chroms}
        self._neg_dict = {key: np.array(dicts["-"][key], dtype="int")
                          for key in self._chroms}
        for val in self._pos_dict.values():
            val.sort()
        for val in self._neg_dict.values():
            val.sort()
        self.total = sum(val.size for val in self._pos_dict.values())
        self.total += sum(val.size for val in self._neg_dict.values())

    def get_chr_names(self):
        return self._chroms
        return ["chr%s" % i for i in range(1, 10)]

    def get_locations_by_chr(self, chrom):
        return self._pos_dict[chrom], self._neg_dict[chrom]

        n_peaks = 10000
        a = np.array([1, 2, 3, 4, 5], dtype="int")
        pos_tags = np.empty(5*n_peaks, dtype="int")
        for i in range(n_peaks):
            pos_tags[(i*5):((i+1)*5)] = i*1000+a
        return pos_tags, pos_tags+50

    @classmethod
    def from_bed_file(cls, file_name):
        dicts = {"+": defaultdict(list),
                 "-": defaultdict(list)}
        with open(file_name) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                chr_name = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[5]
                pos = start if strand == "+" else end
                dicts[strand][chr_name].append(pos)

        return cls(dicts)


class Opt:
    def __init__(self):
        self.umfold = 50
        self.lmfold = 5
        self.bw = 300
        self.info = logging.info
        self.debug = logging.debug
        self.warn = logging.warning
        self.gsize = 300000000


class PeakModel:
    """Peak Model class.
    """

    def __init__(self, opt, treatment, max_pairnum=500):
        self.treatment = treatment
        self.gz = opt.gsize
        self.umfold = opt.umfold
        self.lmfold = opt.lmfold
        self.tag_expansion_size = 10
        self.bw = opt.bw
        self.info = opt.info
        self.debug = opt.debug
        self.warn = opt.warn
        self.error = opt.warn
        self.max_pairnum = max_pairnum
        self.build()

    def build(self):
        """Build the model.
        prepare self.d, self.scan_window, self.plus_line,
        self.minus_line and self.shifted_line to use.
        """
        self.peaksize = 2*self.bw
        self.min_tags = int(round(
            float(self.treatment.total)*self.lmfold*self.peaksize/self.gz/2))
        # mininum unique hits on single strand
        self.max_tags = int(round(
            float(self.treatment.total)*self.umfold*self.peaksize/self.gz/2))
        # maximum unique hits on single strand
        # use treatment data to build model
        self.info("#2 looking for paired plus/minus strand peaks...")
        paired_peakpos = self.__paired_peaks()
        # select up to 1000 pairs of peaks to build model
        num_paired_peakpos = 0
        num_paired_peakpos_remained = self.max_pairnum
        num_paired_peakpos_picked = 0
        # select only num_paired_peakpos_remained pairs.
        for c in paired_peakpos.keys():
            num_paired_peakpos += len(paired_peakpos[c])
        # TL: Now I want to use everything

        num_paired_peakpos_picked = num_paired_peakpos

        self.info("#2 number of paired peaks: %d" % (num_paired_peakpos))
        if num_paired_peakpos < 100:
            logging.warning("Too few paired peaks (%d) to get good shift estimate" % num_paired_peakpos)
            raise NotEnoughPairsException("No enough pairs to build model")
        elif num_paired_peakpos < self.max_pairnum:
            logging.warning("Few paired peaks (%d) than %d! Model may not be build well! Lower your MFOLD parameter may erase this warning. Now I will use %d pairs to build model!" % (num_paired_peakpos, self.max_pairnum,num_paired_peakpos_picked))
        self.__paired_peak_model(paired_peakpos)

    def __paired_peak_model(self, paired_peakpos):
        """Use paired peak positions and treatment
        tag positions to build the model.
        Modify self.(d, model_shift size and
        scan_window size. and extra, plus_line,
        minus_line and shifted_line for plotting).
        """
        window_size = 1+2*self.peaksize+self.tag_expansion_size
        self.plus_line = np.zeros(window_size, dtype="int32")
        self.minus_line = np.zeros(window_size, dtype="int32")
        plus_start = np.zeros(window_size, dtype="int32")
        plus_end = np.zeros(window_size, dtype="int32")
        minus_start = np.zeros(window_size, dtype="int32")
        minus_end = np.zeros(window_size, dtype="int32")
        self.info("start model_add_line...")
        chroms = paired_peakpos.keys()
        for chrom in chroms:
            paired_peakpos_chrom = paired_peakpos[chrom]
            (tags_plus, tags_minus) = self.treatment.get_locations_by_chr(
                chrom)
            self.__model_add_line(paired_peakpos_chrom, tags_plus,
                                  plus_start, plus_end)
            self.__model_add_line(paired_peakpos_chrom, tags_minus,
                                  minus_start, minus_end)

        self.__count(plus_start, plus_end, self.plus_line)
        self.__count(minus_start, minus_end, self.minus_line)

        self.info("start X-correlation...")
        # Now I use cross-correlation to find the best d
        plus_line = self.plus_line
        minus_line = self.minus_line

        # normalize first
        minus_data = (minus_line - minus_line.mean())/(minus_line.std()*len(minus_line))
        plus_data = (plus_line - plus_line.mean())/(plus_line.std()*len(plus_line))

        # cross-correlation
        ycorr = np.correlate(minus_data, plus_data, mode="full")[
            window_size-self.peaksize:window_size+self.peaksize]
        xcorr = np.linspace(len(ycorr)//2*-1, len(ycorr)//2, num=len(ycorr))

        # smooth correlation values to get rid of local maximums from small fluctuations. 
        ycorr = smooth(ycorr, window="flat")

        # all local maximums could be alternative ds.
        i_l_max = np.r_[
            False,
            ycorr[1:] > ycorr[:-1]] & np.r_[ycorr[:-1] > ycorr[1:],
                                            False]
        tmp_cor_alternative_d = ycorr[i_l_max]
        tmp_alternative_d = xcorr[i_l_max]
        cor_alternative_d = tmp_cor_alternative_d[tmp_alternative_d > 0]
        self.alternative_d = [int(val) for val in
                              tmp_alternative_d[tmp_alternative_d > 0]]
        # best cross-correlation point
        self.d = xcorr[np.where(ycorr == max(cor_alternative_d))[0][0]]
        # get rid of the last local maximum if it's at the right end of curve.
        assert len(self.alternative_d) > 0, "No proper d can be found! Tweak --mfold?"

        self.ycorr = ycorr
        self.xcorr = xcorr

        self.scan_window = max(self.d, self.tag_expansion_size)*2
        return True

    def __model_add_line(self, pos1, pos2, start, end):
        """Project each pos in pos2 which is included in
        [pos1-self.peaksize, pos1+self.peaksize] to the line.
        pos1: paired centers -- array.array
        pos2: tags of certain strand -- a numpy.array object
        line: numpy array object where we pileup tags
        """
        psize_adjusted1 = self.peaksize + self.tag_expansion_size // 2

        i1 = 0                  # index for pos1
        i2 = 0                  # index for pos2
        i2_prev = 0             # index for pos2 in previous pos1
        i1_max = len(pos1)
        i2_max = pos2.shape[0]
        flag_find_overlap = False

        max_index = start.shape[0] - 1

        while i1 < i1_max and i2 < i2_max:
            p1 = pos1[i1]
            p2 = pos2[i2]
            if p1-psize_adjusted1 > p2:
                i2 += 1
            elif p1 + psize_adjusted1 < p2:
                i1 += 1
                i2 = i2_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    i2_prev = i2
                s = max(p2-p1+self.peaksize, 0)
                start[s] += 1
                e = min(p2+self.tag_expansion_size-p1+self.peaksize,
                        max_index)
                end[e] -= 1
                i2 += 1
        return

    def __model_add_line_vanilla(self, pos1, pos2, start, end):
        """Project each pos in pos2 which is included in
        [pos1-self.peaksize, pos1+self.peaksize] to the line.
        pos1: paired centers -- array.array
        pos2: tags of certain strand -- a numpy.array object
        line: numpy array object where we pileup tags
        """
        i1 = 0                  # index for pos1
        i2 = 0                  # index for pos2
        i2_prev = 0             # index for pos2 in previous pos1
        i1_max = len(pos1)
        i2_max = pos2.shape[0]
        flag_find_overlap = False

        max_index = start.shape[0] - 1
        psize_adjusted1 = self.peaksize + self.tag_expansion_size // 2

        while i1 < i1_max and i2 < i2_max:
            p1 = pos1[i1]
            p2 = pos2[i2]
            if p1-psize_adjusted1 > p2:
                i2 += 1
            elif p1+psize_adjusted1 < p2:
                i1 += 1
                i2 = i2_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    i2_prev = i2
                s = max(p2-self.tag_expansion_size//2-p1+psize_adjusted1, 0)
                start[s] += 1
                e = min(p2+self.tag_expansion_size//2-p1+psize_adjusted1,
                        max_index)
                end[e] -= 1
                i2 += 1
        return

    def __count(self, start, end, line):
        pileup = 0
        for i in range(line.shape[0]):
            pileup += start[i] + end[i]
            line[i] = pileup
        return

    def __paired_peaks(self):
        """Call paired peaks from fwtrackI object.
        Return paired peaks center positions.
        """
        chrs = self.treatment.get_chr_names()
        chrs.sort()
        paired_peaks_pos = {}
        for i in range(len(chrs)):
            chrom = chrs[i]
            plus_tags, minus_tags = self.treatment.get_locations_by_chr(chrom)
            plus_peaksinfo = self.__naive_find_peaks(plus_tags, 1)
            minus_peaksinfo = self.__naive_find_peaks(minus_tags, 0)
            if not plus_peaksinfo or not minus_peaksinfo:
                continue
            paired_peaks_pos[chrom] = self.__find_pair_center(
                plus_peaksinfo, minus_peaksinfo)
        return paired_peaks_pos

    def __find_pair_center(self, pluspeaks, minuspeaks):
        ip = 0                  # index for plus peaks
        im = 0                  # index for minus peaks
        im_prev = 0             # index for minus peaks in previous plus peak
        pair_centers = []
        ip_max = len(pluspeaks)
        im_max = len(minuspeaks)
        flag_find_overlap = False
        while ip < ip_max and im < im_max:
            (pp, pn) = pluspeaks[ip]  # for (peakposition, tagnumber in peak)
            (mp, mn) = minuspeaks[im]
            if pp-self.peaksize > mp:  # move minus
                im += 1
            elif pp+self.peaksize < mp:  # move plus
                ip += 1
                im = im_prev    # search minus peaks from previous index
                flag_find_overlap = False
            else:               # overlap!
                if not flag_find_overlap:
                    flag_find_overlap = True
                    im_prev = im  # only the first index is recorded
                # number tags in plus and minus peak region are comparable...
                if float(pn)/mn < 2 and float(pn)/mn > 0.5:
                    if pp < mp:
                        pair_centers.append((pp+mp)//2)
                im += 1
        return pair_centers

    def __naive_find_peaks(self, taglist, plus_strand=1):
        """Naively call peaks based on tags counting.
        if plus_strand == 0, call peak on minus strand.
        Return peak positions and the tag number in peak
        region by a tuple list [(pos,num)].
        """
        peak_info = []    # store peak pos in every peak region and
                          # unique tag number in every peak region
        if taglist.shape[0] < 2:
            return peak_info
        pos = taglist[0]
        # current_tag_list = [pos]
        cur_start = 0
        for i in range(1, len(taglist)):
            pos = taglist[i]
            # call peak in current_tag_list when the region is long enough
            if (pos - taglist[cur_start] + 1) > self.peaksize:
                # a peak will be called if tag number is ge min tags.
                if self.max_tags >= i-cur_start >= self.min_tags:
                    peak_info.append((self.__naive_peak_pos(
                        taglist[cur_start:i], plus_strand),
                                      i-cur_start))
                cur_start = i

            # current_tag_list.append(pos)   # add pos while 1. no
                                             # need to call peak;
                                             # 2. current_tag_list is []
        return peak_info

    def __get_horizon_line_sparse(self, pos_list, start, peak_length):
        positions = np.empty(pos_list.size*2)
        positions[:pos_list.size] = np.maximum(
            pos_list-start-self.tag_expansion_size//2,
            0)
        positions[pos_list.size:] = np.minimum(
            pos_list-start+self.tag_expansion_size//2,
            peak_length)
        args = np.argsort(positions)
        codes = 1-2*(args >= pos_list.size)

    def __get_horizon_line(self, pos_list):
        # pos_list = np.array(pos_list, dtype="int")
        peak_length = pos_list[-1]+1-pos_list[0] + self.tag_expansion_size
        # leftmost position of project line
        start = pos_list[0] - self.tag_expansion_size//2
        ss = pos_list-start-self.tag_expansion_size//2
        es = pos_list-start+self.tag_expansion_size//2

        # the line for tags to be projected
        horizon_line = np.zeros(peak_length, dtype="int32")
        horizon_line[ss] += 1
        horizon_line[es] -= 1
        horizon_line = np.cumsum(horizon_line)
        return horizon_line

    def __naive_peak_pos(self, pos_list, plus_strand):
        """Naively calculate the position of peak.
        plus_strand: 1, plus; 0, minus
        return the highest peak summit position.
        """
        peak_length = pos_list[-1]+1-pos_list[0] + self.tag_expansion_size
        start = pos_list[0] - self.tag_expansion_size//2
        horizon_line = self.__get_horizon_line(pos_list)
        top_pos = self.__find_top_pos(horizon_line, peak_length)
        return (top_pos[len(top_pos)//2]+start)


    def __find_top_pos(self, horizon_line, peak_length):
        m = np.max(horizon_line)
        return np.where(horizon_line == m)[0]



# smooth function from SciPy cookbook: http://www.scipy.org/Cookbook/SignalSmooth
def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the beginning and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
        
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[(window_len//2):-(window_len//2)]

