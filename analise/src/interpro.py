# https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/P04637
# https://github.com/ProteinsWebTeam/interpro7-api/tree/master/docs
import json
import os

import requests

from src.structure_data import RegionData


class InterproData:
    def __init__(self, database, output_dir):
        self.database = database
        self.output_dir = output_dir
        self.parsed_data = {}

    def prepare_interpro_data_to_download(self):
        # to tak nie działa
        uniprot_ids = []
        for virus, data in self.database.database.items():
            for lcr in data:
                uniprot_ids.append(lcr.uniprot_id)
        uniprot_ids = list(set(uniprot_ids))
        for uniprot in uniprot_ids:
            url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/" + uniprot
            print(url)
            r = requests.get(url)
            with open(self.output_dir + uniprot + ".json", "w") as f:
                f.write(r.text)

    def parse_files(self):
        files = os.listdir(self.output_dir)
        for file_name in files:
            with open(self.output_dir + file_name) as f:
                print(self.output_dir + file_name, f.name)
                file_data = f.read()
                if not file_data:
                    continue
                soup = json.loads(file_data)
                self.parsed_data[file_name.split(".")[0]] = []
                for result in soup['results']:
                    meta = result["metadata"]
                    proteins = result["proteins"][0]
                    acc = meta["accession"]
                    name_ip = meta["name"]
                    source = meta["source_database"]
                    type_data = meta["type"]
                    member_database = meta["member_databases"]
                    for database, data_database in member_database.items():
                        if data_database is not None:
                            for database_id, name in data_database.items():
                                if isinstance(name, list):
                                    pass
                                else:
                                    uniprot = proteins["accession"].upper()
                                    regions = proteins["entry_protein_locations"]
                                    for region in regions:
                                        for fragment in region["fragments"]:
                                            start = fragment["start"]
                                            end = fragment["end"]
                                            desc = ",".join([acc, name_ip, type_data, database_id, name])
                                            self.parsed_data[uniprot].append(
                                                RegionData(uniprot, start, end, desc, source + "/" + database))
