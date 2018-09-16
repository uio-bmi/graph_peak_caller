import sys
from graph_peak_caller.analysis.peakscomparer import AnalysisResults
import logging
import numpy as np


class HtmlReportGenerator:
    def __init__(self, transcription_factors):
        self.tfs = transcription_factors
        self.html = ""
        self.latex_table = ""
        #self._create_report_table()
        self._create_simple_report_table()
        self.create_bar_plots()
        self._create_report_figures()

    def _write_simple_table_row(self, tf, analysis_result):
        assert analysis_result.tot_peaks1 == analysis_result.tot_peaks2, "Only use simple report when number of peaks are the same"
        assert analysis_result.peaks1_in_peaks2 == analysis_result.peaks2_in_peaks1, "Only use simple report when number of peaks are the same"
        self.html += "<tr>"
        self.html += "<td>%s</td>" % tf
        self.html += "<td>%d</td>" % (analysis_result.tot_peaks1)
        self.html += "<td>%d</td>" % (analysis_result.peaks1_in_peaks2)
        self.html += "<td>%d</td>" % (analysis_result.peaks1_not_in_peaks2)
        self.html += "<td>%d (%.3f %%)</td>" % (analysis_result.peaks1_in_peaks2_matching_motif, 100 * analysis_result.peaks1_in_peaks2_matching_motif / analysis_result.peaks1_in_peaks2)
        self.html += "<td>%d (%.3f %%)</td>" % (analysis_result.peaks1_not_in_peaks2_matching_motif, 100 * analysis_result.peaks1_not_in_peaks2_matching_motif / analysis_result.peaks1_not_in_peaks2)
        self.html += "<td>%d (%.3f %%)</td>" % (analysis_result.peaks2_in_peaks1_matching_motif, 100 * analysis_result.peaks2_in_peaks1_matching_motif / analysis_result.peaks2_in_peaks1)
        self.html += "<td>%d (%.3f %%)</td>" % (analysis_result.peaks2_not_in_peaks1_matching_motif, 100 * analysis_result.peaks2_not_in_peaks1_matching_motif / analysis_result.peaks2_not_in_peaks1)
        self.html += "<td>%.3f</td>" % np.mean(analysis_result.peaks1_in_peaks2_bp_not_on_linear)
        self.html += "<td>%.3f</td>" % np.mean(analysis_result.peaks1_not_in_peaks2_bp_not_on_linear)
        self.html += "</tr>"
                
        #   \multicolumn{1}{r|}{SOC1} &15553 & 1664/14314  (11.6\%) & 117/1239    (9.44\%) &16407 &1681/14314   (11.7\%)  &142/2093  (6.78\%) \\
        self.latex_table += "\multicolumn{1}{r|}{" + tf.split("_")[-1] + "} & "
        self.latex_table += "%d & " % analysis_result.tot_peaks1
        self.latex_table += "%d & " % analysis_result.peaks1_in_peaks2
        self.latex_table += "%d & " % analysis_result.peaks1_not_in_peaks2
        self.latex_table += "\multicolumn{1}{|l}{%d} & %.2f\%% & " % (analysis_result.peaks1_in_peaks2_matching_motif, 100 * analysis_result.peaks1_in_peaks2_matching_motif / analysis_result.peaks1_in_peaks2)
        self.latex_table += "%d & %.2f\%% & " % (analysis_result.peaks1_not_in_peaks2_matching_motif, 100 * analysis_result.peaks1_not_in_peaks2_matching_motif / analysis_result.peaks1_not_in_peaks2)
        self.latex_table += "\multicolumn{1}{|l}{%d} & %.2f\%% & " % (analysis_result.peaks2_in_peaks1_matching_motif, 100 * analysis_result.peaks2_in_peaks1_matching_motif / analysis_result.peaks2_in_peaks1)
        self.latex_table += "%d & %.2f\%%" % (analysis_result.peaks2_not_in_peaks1_matching_motif, 100 * analysis_result.peaks2_not_in_peaks1_matching_motif / analysis_result.peaks2_not_in_peaks1)
        self.latex_table += "\\\ \n" 
        
    
    def _write_table_row(self, tf, analysis_result):
        self.html += """
        <tr>
            <td>%s</td>
            <td>%d</td>
            <td>%d (%d) = %.3f </td>
            <td>%d (%d) = %.3f</td>
            <td>%d</td>
            <td>%d (%d) = %.3f</td>
            <td>%d (%d) = %.3f</td>
            <td>%.2f</td>
            <td>%.2f</td>
            <td>%.4f</td>
            <td>%.4f</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%d</td>
            <td>%d</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.3f</td>
            <td>%.3f</td>
        </tr>
        """ % (tf,
               analysis_result.tot_peaks1,
               analysis_result.peaks1_in_peaks2,
               analysis_result.peaks1_in_peaks2_matching_motif,
               100 * analysis_result.peaks1_in_peaks2_matching_motif / analysis_result.peaks1_in_peaks2,
               analysis_result.peaks1_not_in_peaks2,
               analysis_result.peaks1_not_in_peaks2_matching_motif,
               100 * analysis_result.peaks1_not_in_peaks2_matching_motif / analysis_result.peaks1_not_in_peaks2,
               analysis_result.tot_peaks2,
               analysis_result.peaks2_in_peaks1,
               analysis_result.peaks2_in_peaks1_matching_motif,
               100 * analysis_result.peaks2_in_peaks1_matching_motif / analysis_result.peaks2_in_peaks1,
               analysis_result.peaks2_not_in_peaks1,
               analysis_result.peaks2_not_in_peaks1_matching_motif,
               100 * analysis_result.peaks2_not_in_peaks1_matching_motif / analysis_result.peaks2_not_in_peaks1,
               np.mean(analysis_result.peaks1_in_peaks2_bp_not_on_linear),
               np.mean(analysis_result.peaks1_not_in_peaks2_bp_not_on_linear),
               analysis_result.peaks1_total_nodes / analysis_result.peaks1_total_basepairs,
               analysis_result.peaks2_total_nodes / analysis_result.peaks2_total_basepairs,
               analysis_result.peaks1_total_basepairs / analysis_result.tot_peaks1,
               analysis_result.peaks2_total_basepairs / analysis_result.tot_peaks2,
               analysis_result.peaks1_total_basepairs,
               analysis_result.peaks2_total_basepairs,
               analysis_result.peaks1_unique_total_basepairs,
               analysis_result.peaks2_unique_total_basepairs,
               analysis_result.peaks1_unique_total_nodes / analysis_result.peaks1_unique_total_basepairs,
               analysis_result.peaks2_unique_total_nodes / analysis_result.peaks2_unique_total_basepairs,
               np.mean(analysis_result.peaks1_unique_scores),
               np.mean(analysis_result.peaks2_unique_scores),
               )

    def _create_simple_report_table(self):                                                                                                                                                                                                    
        self.html += """
        <table class='table' style='margin-top: 40px;'>
            <h3>Overview of peaks found.</h3>
            <thead>
                <tr>
                    <th></th>
                    <th colspan='3'>Peaks found</th>
                    <th colspan='2'>Graph Peak Caller motif enrichment</th>
                    <th colspan='2'>MACS motif enrichment</th>
                    <th colspan='2'>BP not on ref</th>
                </tr>
            </theads>
            <tr>
                <th>TF</th>
                <th># Found in total</th>
                <th># Found by both</th>
                <th># Unique</th>
                <th># Among shared peaks</th>
                <th># Among unique peaks</th>
                <th># Among shared peaks</th>
                <th># Among unique peaks</th>
                <th>GPC</th>
                <th>MACS</th>
            </tr>
        """
        summed_results = AnalysisResults()
        for tf in self.tfs:
            results = AnalysisResults.from_file("figures_tables/" + tf + ".pickle")
            self._write_simple_table_row(tf, results)
            summed_results += results

        self._write_simple_table_row("SUM", summed_results)

        self.html += "</table>"

        self.html += "<pre>" + self.latex_table + "</pre>"

    def _create_report_table(self):
        self.html += """
        <table class='table' style='margin-top: 40px;'>
            <h3>Overview of peaks found.</h3>
            <thead>
                <tr>
                    <th></th>
                    <th colspan='3'>Graph Peak Caller</th>
                    <th colspan='3'>Macs2</th>
                    <th colspan='4'>Graph Peak Caller motif enrichment</th>
                </tr>
            </theads>
            <tr>
                <th>TF</th>
                <th># Peaks found</th>
                <th># Peaks also found by MACS2</th>
                <th># Peaks NOT found by MACS2</th>
                <th># Peaks found</th>
                <th># Peaks also found by Graph Peak Caller</th>
                <th># Peaks NOT found by Graph Peak Caller</th>
                <th>Average number of base pairs of GPC peaks also found by MACS2 that are part of linear reference genome</th>
                <th>Average number of base pairs of GPC peaks NOT found by MACS2 that are part of linear reference genome</th>
                <th>Nodes per base pair graph peak caller</th>
                <th>Nodes per base pair MACS2</th>
                <th>Average size gpc</th>
                <th>Average size macs2</th>
                <th>Total base pairs gpc</th>
                <th>Total base pairs macs</th>
                <th>Unique gpc peaks total base pairs</th>
                <th>Unique macs peaks total base pairs</th>
                <th>Unique gpc peaks nodes per bp</th>
                <th>Unique macs peaks nodes per bp</th>
                <th>Avg score unique gpc</th>
                <th>Avg score unique macs</th>
            </tr>
        """
        summed_results = AnalysisResults()
        for tf in self.tfs:
            results = AnalysisResults.from_file("figures_tables/" + tf + ".pickle")
            self._write_table_row(tf, results)
            summed_results += results

        self._write_table_row("SUM", summed_results)

        self.html += "</table>"

    def _create_report_figures(self):
        self.html += "<h3>Motifenrichment plots</h3>"
        for tf in self.tfs:
            self.html += "<p>All peaks.</p>" \
                         "<img style='width: 700px; height: auto; padding: 50px;' src='" + tf + ".png'/>"
            self.html += "<p>Unique peaks:</p>" \
                         "<img style='width: 700px; height: auto; padding: 50px;' src='" + tf + "_unique_peaks.png'/>"

    def _html_start(self):
        import datetime
        return """
        <!doctype html>
        <html>
        <head>
            <link rel="stylesheet" href="bootstrap.min.css">
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
            <script src="https://code.highcharts.com/highcharts.js"></script>
            <script src="https://code.highcharts.com/modules/exporting.js"></script>
            <script src="https://code.highcharts.com/modules/export-data.js"></script>

            <title>Graph Peak Caller - Jenkins results</title>
        </head>
        <body>
        <div class="container" style='margin-top: 50px; max-width: 1600px;'>


        <h1>Graph Peak Caller - experiment results</h1>
        <h2>Report generated %s</h2>
        <div id='bar-plot'></div>
        """ % datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")

    def _html_end(self):
        return "</div></body></html>"

    def write_report(self, file_name):
        f = open(file_name, "w")
        final_html = self._html_start() + self.html + self._html_end()
        f.write(final_html)
        f.close()

    def create_bar_plots(self):

        category_names = "'" + "', '".join([tf.replace("ARABIDOPSIS_", "") for tf in self.tfs]) + "'"

        self.html += """
            <script>
            Highcharts.chart('bar-plot', {

                chart: {
                    type: 'column'
                },

                title: {
                    text: ''
                },

                xAxis: {
                    categories: [%s]
                },

                yAxis: {
                    allowDecimals: false,
                    min: 0,
                    title: {
                        text: 'Number of peaks'
                    }
                },

                tooltip: {
                    formatter: function () {
                        return '<b>' + this.x + '</b><br/>' +
                            this.series.name + ': ' + this.y + '<br/>' +
                            'Total: ' + this.point.stackTotal;
                    }
                },

                plotOptions: {
                    column: {
                        stacking: 'percent'
                    }
                },
        """ % (category_names)

        total_gpc = []
        total_gpc_with_motif = []
        total_macs = []
        total_macs_with_motif = []
        for tf in self.tfs:
            results = AnalysisResults.from_file("figures_tables/" + tf + ".pickle")
            total_gpc.append(results.peaks1_not_in_peaks2)
            total_macs.append(results.peaks2_not_in_peaks1)
            total_gpc_with_motif.append(results.peaks1_not_in_peaks2_matching_motif)
            total_macs_with_motif.append(results.peaks2_not_in_peaks1_matching_motif)



        self.html += """series: [

                {
                    name: 'MACS2 unique peaks',
                    data: [""" + ','.join(str(t) for t in total_macs) + """],
                    stack: 'male'
                }, {
                    name: 'MACS2 unique peaks with motif hit',
                    data: [""" + ','.join(str(t) for t in total_macs_with_motif) + """],
                    stack: 'male'
                }, {
                    name: 'Graph Peak Caller unique peaks',
                    data: [""" + ','.join(str(t) for t in total_gpc) + """],
                    stack: 'female'
                }, {
                    name: 'Graph Peak Caller unique peaks with motif hit',
                    data: [""" + ','.join(str(t) for t in total_gpc_with_motif) + """],
                    stack: 'female'
                }
        ]
            });
            </script>
            """

if __name__ == "__main__":
    transcription_factors = sys.argv[1].split(",")
    report_file_name = sys.argv[2]

    generator = HtmlReportGenerator(transcription_factors)
    generator.write_report(report_file_name)
