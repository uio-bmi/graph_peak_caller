# -*- coding: latin-1 -*-
from urllib import request
import json
import sys
import requests, json
import gzip
import os.path
import logging
import shutil
import wget

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}
BASE_URL = "https://www.encodeproject.org/"


class EncodeWrapper(object):
    def __init__(self):
        self.BASE_URL = BASE_URL

    def _call_api(self, url):
        url = self.BASE_URL + url
        response = requests.get(url, headers=HEADERS)
        if response.status_code != 200:
            raise Exception("Error while calling API %s . Error: %s"
                            % (url, ""))

        return response.json()

    def download_file(self, data_url, out_file):

        print("Downloading file %s" % data_url)

        if os.path.isfile(out_file):
            logging.info("%s has been downloaded before, not downloading again"
                         % data_url)
            return

        logging.info("Downloading file %s to %s" % (data_url, out_file))
        wget.download(BASE_URL + data_url)


def get_experiment_raw_data_url(accession_number, replicate_number=1):
    """
    Downloads fastq file for given accession ID and replicate nuber
    """
    encode = EncodeWrapper()
    data = encode._call_api("experiments/%s" % accession_number)

    file_url = None

    for file in data["files"]:
        if file["output_category"] == "raw data":
            #print("File %s with replicate number %s" % (file["href"], file["replicate"]["biological_replicate_number"]))
            #print(file["file_type"])
            #print(replicate_number)
            if file["file_type"] == "fastq" and \
                    str(file["replicate"]["biological_replicate_number"]) == str(replicate_number):
                file_url = file["href"]

    assert file_url is not None
    print(BASE_URL + file_url)
    #file_name = "raw.fastq.gz"
    #out_file = "tmp/%s.gz" % file_name
    #encode.download_file(file_url, out_file)


if __name__ == "__main__":

    encode_id = sys.argv[1]
    replicate_number = sys.argv[2]
    get_experiment_raw_data_url(encode_id, replicate_number)

