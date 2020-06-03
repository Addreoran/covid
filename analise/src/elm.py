# https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchRPSBWebAPI_GETorPOST
# https://www.ncbi.nlm.nih.gov/protein/P98073
# https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?queries=P98073&useid1=true&tdata=hits&db=cdd,pfam,smart,cog,tigrfam,kog&smode=live&filter=false&qdefl=true&cddefl=true
from src.structure_data import RegionData


class ELMInstance:
    def __init__(self, uniprot_idies, begin, end, accession, ELMIdentifier):
        self.uniprot_idies = uniprot_idies
        self.begin = begin
        self.end = end
        self.accession = accession
        self.ELMIdentifier = ELMIdentifier


class ELMDatabase:
    def __init__(self, database, output_dir):
        self.database = database
        self.output_dir = output_dir
        self.parsed_data = {}

    def prepare_elm_data(self, path):
        elm_database = []
        with open(path) as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    if not line.startswith('"Accession'):
                        line = [i.replace('"', "") for i in line.split("\t")]
                        elm_database.append(ELMInstance(line[5].split(), line[6], line[7], line[0], line[1]))
        self.select_data(elm_database)

    def select_data(self, elm_database):
        uniprot_ids = []
        for virus, data in self.database.database.items():
            for lcr in data:
                uniprot_ids.append(lcr.uniprot_id)
        uniprot_ids = list(set(uniprot_ids))
        for uniprot in uniprot_ids:
            self.parsed_data[uniprot] = []
            for elm_ins in elm_database:
                if uniprot in elm_ins.uniprot_idies:
                    self.parsed_data[uniprot].append(RegionData(uniprot, elm_ins.begin, elm_ins.end,
                                                                elm_ins.accession + "," + elm_ins.ELMIdentifier, "elm"))

