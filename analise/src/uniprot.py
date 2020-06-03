import os

import requests
from bs4 import BeautifulSoup

from src.structure_data import RegionData


class UniprotData:
    def __init__(self, database, output_dir):
        self.database = database
        self.output_dir = output_dir
        self.parsed_data = {}

    def prepare_uniprot_data_to_download(self):
        uniprot_ids = []
        for virus, data in self.database.database.items():
            for lcr in data:
                uniprot_ids.append(lcr.uniprot_id)
        uniprot_ids = list(set(uniprot_ids))
        for uniprot in uniprot_ids:
            url = "https://www.uniprot.org/uniprot/" + uniprot + ".xml"
            print(url)
            r = requests.get(url)
            with open(self.output_dir + uniprot + ".xml", "w") as f:
                f.write(r.text)

    def parse_files(self):
        #     def __init__(self, begin, id, end, desc, database):
        # tmp = RegionData()
        files = os.listdir(self.output_dir)
        for file_name in files:
            with open(self.output_dir + file_name) as f:
                soup = BeautifulSoup(f.read(), 'html.parser')
                features = soup.find_all("feature")
                uniprot = file_name.split(".")[0]
                self.parsed_data[uniprot] = []
                for feature in features:
                    begin = feature.find("begin")
                    end = feature.find("end")
                    position = None
                    if begin is None:
                        position = feature.find("position")
                        begin = position
                        end = position
                    try:
                        descr = feature['description']
                    except:
                        descr = None
                    try:
                        id_desc = feature["id"]
                    except:
                        id_desc = None
                    try:
                        type = feature['type']
                    except:
                        type = None
                    try:
                        begin = begin['position']
                        end = end['position']
                        desc = ",".join([str(descr), str(id_desc), str(type)])
                        self.parsed_data[uniprot].append(
                            RegionData(uniprot, begin, end, desc, "uniprot"))
                    except:
                        print(file_name)
