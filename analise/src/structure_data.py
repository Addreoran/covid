class LCRData:
    def __init__(self, header, virus):
        self.header = header
        self.virus = virus
        self.uniprot_id = None
        self.begin = None
        self.end = None
        # {region:{lcr_cover:10, region_cover:10}}
        self.cover = None

    def set_uniprot_id(self):
        self.uniprot_id = self.header.split("|")[1]
        self.begin = int(self.header.split(":")[-1].split(",")[0].replace("begin=", ""))
        self.end = int(self.header.split(":")[-1].split(",")[1].split("Length=")[0].replace("end=", ""))
        print(self.begin, self.end)


class RegionData:
    def __init__(self, id, begin, end, desc, database):
        self.begin = int(begin)
        self.end = int(end)
        self.id = id
        self.desc = desc
        self.database = database
