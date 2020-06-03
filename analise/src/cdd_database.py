# https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchRPSBWebAPI_GETorPOST
# https://www.ncbi.nlm.nih.gov/protein/P98073
# https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?queries=P98073&useid1=true&tdata=hits&db=cdd,pfam,smart,cog,tigrfam,kog&smode=live&filter=false&qdefl=true&cddefl=true
import os
import time

import requests

from src.structure_data import RegionData


class CDDData:
    def __init__(self, database, output_dir):
        self.database = database
        self.output_dir = output_dir
        self.parsed_data = {}

    def prepare_cdd_data_to_download(self):
        uniprot_ids = []
        for virus, data in self.database.database.items():
            for lcr in data:
                uniprot_ids.append(lcr.uniprot_id)
        uniprot_ids = list(set(uniprot_ids))
        urls_vs = {}
        for uniprot in uniprot_ids:
            url = f"https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?queries={uniprot}&useid1=true&tdata=hits&db=cdd,pfam,smart,cog,tigrfam,kog&smode=live&filter=false&qdefl=true&cddefl=true"
            print(url)
            r = requests.get(url)
            for i in r.text.splitlines():
                if "#cdsid" in i:
                    cdsid = i.split()[-1]
                    print("nr", i)
            url2 = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?cdsid=" + cdsid
            urls_vs[url2] = uniprot
        time.sleep(20)

        for url2, uniprot in urls_vs.items():
            r = requests.get(url2)
            print(url2)
            with open(self.output_dir + uniprot + ".csv", "w") as f:
                f.write(r.text)

    def parse_files(self):
        files = os.listdir(self.output_dir)
        for file_name in files:
            test = False
            with open(self.output_dir + file_name) as f:
                uniprot = file_name.split(".")[0]
                self.parsed_data[uniprot] = []
                for line in f.readlines():
                    # print(line)
                    line = line.strip()
                    if line.startswith("Query"):
                        test = True
                    elif test:
                        line = line.split("\t")
                        desc = ",".join([line[7], line[8], line[-1]])
                        print((uniprot, line[3], line[4], desc, "cd-search"))
                        self.parsed_data[uniprot].append(
                            RegionData(uniprot, line[3], line[4], desc, "cd-search"))

                    # soup = BeautifulSoup(f.read(), 'html.parser')
                    # features = soup.find_all("match")
                    # uniprot = file_name.split(".")[0]
                    # self.parsed_data[uniprot] = []
                    # for feature in features:
                    #     begin = feature.find("start").text
                    #     end = feature.find("stop").text
                    #     signature_ac = feature.find("signature_ac").text
                    #     signature_id = feature.find("signature_id").text
                    #     print(begin, end, signature_ac, signature_id)
                    #     desc = ",".join([str(signature_id), str(signature_ac)])
                    #     self.parsed_data[uniprot].append(
                    #         RegionData(uniprot, begin, end, desc, "prosite"))
