# https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi?seq=P98073 do znalezienia motywów
# https://prosite.expasy.org/scanprosite/scanprosite_doc.html
# https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi?seq=P98073&output=xml
import os

import requests
from bs4 import BeautifulSoup

from src.structure_data import RegionData


class PrositeData:
    def __init__(self, database, output_dir):
        self.database = database
        self.output_dir = output_dir
        self.parsed_data = {}

    def prepare_prosite_data_to_download(self):
        downloaded = os.listdir(self.output_dir)
        downloaded = [i.split(".")[0] for i in downloaded]
        uniprot_ids = []
        for virus, data in self.database.database.items():
            for lcr in data:
                uniprot_ids.append(lcr.uniprot_id)
        uniprot_ids = list(set(uniprot_ids))
        for uniprot in uniprot_ids:
            url = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi?seq=" + uniprot + "&output=xml"
            print(url)
            if uniprot not in downloaded:
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
                features = soup.find_all("match")
                uniprot = file_name.split(".")[0]
                self.parsed_data[uniprot] = []
                for feature in features:
                    begin = feature.find("start").text
                    end = feature.find("stop").text
                    signature_ac = feature.find("signature_ac").text
                    signature_id = feature.find("signature_id").text
                    print(begin, end, signature_ac, signature_id)
                    desc = ",".join([str(signature_id), str(signature_ac)])
                    self.parsed_data[uniprot].append(
                        RegionData(uniprot, begin, end, desc, "prosite"))
