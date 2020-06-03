from src.cdd_database import CDDData
from src.elm import ELMDatabase
from src.interpro import InterproData
from src.prosite import PrositeData
from src.structure_data import LCRData
from src.uniprot import UniprotData


class LoadDatabase:
    def __init__(self, input_file, output_dir):
        self.input = input_file
        self.output_dir = output_dir
        self.database = {}
        self.regions_info = {}

    def read_input(self):
        with open(self.input) as f:
            for line in f.readlines():
                if not line.startswith("query_header"):
                    line = line.split(";")
                    if line[0] not in self.database:
                        self.database[line[0]] = []
                    tmp = LCRData(line[1], line[0])
                    tmp.set_uniprot_id()
                    self.database[line[0]].append(tmp)

    def read_gbsc_input(self):
        viruses = []
        save = True
        regions = []
        with open(self.input) as f:
            for line in f.readlines():
                # print(line)
                line = line.strip()
                if line.startswith("Human"):
                    print("tu")
                    save = False
                elif not line.startswith("SARS") and save:
                    viruses.append(line.split(":", 1)[1].replace('"', ""))

                elif not save:
                    print(line.split(":", 1), line)
                    regions.append(line.split(":", 1)[1].replace('"', ""))

        print(viruses)
        print(regions)
        for virus in viruses:
            self.database[virus] = []
            for region in regions:
                tmp = LCRData(region, virus)
                tmp.set_uniprot_id()
                self.database[virus].append(tmp)

    def add_region_info(self, data):
        for protein, regions in data.items():
            if protein not in self.regions_info.keys():
                self.regions_info[protein] = regions
            else:
                self.regions_info[protein] += regions
                self.regions_info[protein] = list(set(self.regions_info[protein]))

    def check_covarage_with_lcr(self):
        for virus, lcrs in self.database.items():
            for lcr in lcrs:
                lcr.cover = {}
                for functional_region in self.regions_info[lcr.uniprot_id]:
                    lcr_cover, region_cover = self.check_if_cover(lcr, functional_region)
                    print(lcr.begin, lcr.end, functional_region.begin, functional_region.end, lcr_cover, region_cover)
                    lcr.cover[functional_region] = dict(lcr_cover=lcr_cover, region_cover=region_cover)

    def make_summary(self, file):
        # białko wirusowe | funkcja białka wirusowego | subunit | funkcja białka ludzkiego | ile białek ludzkich ma funkcję | 1-20% | 21-40% | 41-60% | 61-80% | 81-100% | ile jest tych białek ludzkich
        for virus, lcrs in self.database.items():
            functional = {}
            proteins = []
            for lcr in lcrs:
                proteins.append(lcr.uniprot_id)
                for functions, covers in lcr.cover.items():
                    if (functions.desc, functions.database) not in functional.keys():
                        functional[(functions.desc, functions.database)] = {}
                        functional[(functions.desc, functions.database)]["lcr_cover"] = {(0, 10): 0, (10, 20): 0,
                                                                                         (20, 30): 0,
                                                                                         (30, 40): 0, (40, 50): 0,
                                                                                         (50, 60): 0,
                                                                                         (60, 70): 0, (70, 80): 0,
                                                                                         (80, 90): 0,
                                                                                         (90, 100): 0}
                        functional[(functions.desc, functions.database)]["region_cover"] = {(0, 10): 0, (10, 20): 0,
                                                                                            (20, 30): 0,
                                                                                            (30, 40): 0, (40, 50): 0,
                                                                                            (50, 60): 0,
                                                                                            (60, 70): 0, (70, 80): 0,
                                                                                            (80, 90): 0,
                                                                                            (90, 100): 0}
                    for cat, value in covers.items():
                        for i in functional[(functions.desc, functions.database)][cat].keys():
                            if i[0] < value <= i[1]:
                                functional[(functions.desc, functions.database)][cat][(i[0], i[1])] += 1
            with open(file, "a") as f:
                f.write(
                    "virus?function?database?SUm of LCR cover?"
                    "LCR cover: 0-10?LCR cover: 11-20?LCR cover: 21-30?LCR cover: 31-40?LCR cover: 41-50?"
                    "LCR cover: 51-60?LCR cover: 61-70?LCR cover: 71-80?LCR cover: 81-90?LCR cover: 91-100?"
                    "Sum of region cover?Region cover: 0-10?Region cover: 11-20?Region cover: 21-30?Region cover: 31-40?"
                    "Region cover: 41-50?Region cover: 51-60?Region cover: 61-70?Region cover: 71-80?"
                    "Region cover: 81-90?Region cover: 91-100?all LCR regions?number of all proteins\n")

                for data in functional.items():
                    # print(list(data[1]["lcr_cover"]))
                    if sum(list(data[1]["lcr_cover"].values())) > 0 or sum(
                            list(data[1]["region_cover"].values())) > 0:
                        f.write("?".join([str(i) for i in [virus,
                                                           data[0][0],
                                                           data[0][1],
                                                           sum(list(data[1]["lcr_cover"].values())),
                                                           data[1]["lcr_cover"][(0, 10)],
                                                           data[1]["lcr_cover"][(10, 20)],
                                                           data[1]["lcr_cover"][(20, 30)],
                                                           data[1]["lcr_cover"][(30, 40)],
                                                           data[1]["lcr_cover"][(40, 50)],
                                                           data[1]["lcr_cover"][(50, 60)],
                                                           data[1]["lcr_cover"][(60, 70)],
                                                           data[1]["lcr_cover"][(70, 80)],
                                                           data[1]["lcr_cover"][(80, 90)],
                                                           data[1]["lcr_cover"][(90, 100)],
                                                           sum(data[1]["region_cover"].values()),
                                                           data[1]["region_cover"][(0, 10)],
                                                           data[1]["region_cover"][(10, 20)],
                                                           data[1]["region_cover"][(20, 30)],
                                                           data[1]["region_cover"][(30, 40)],
                                                           data[1]["region_cover"][(40, 50)],
                                                           data[1]["region_cover"][(50, 60)],
                                                           data[1]["region_cover"][(60, 70)],
                                                           data[1]["region_cover"][(70, 80)],
                                                           data[1]["region_cover"][(80, 90)],
                                                           data[1]["region_cover"][(90, 100)],
                                                           len(self.database[virus]),
                                                           len(list(set(proteins))), "\n"]]))

                        pass

    def check_if_cover(self, lcr, region):
        lcr_cover, region_cover = 0, 0
        # ----|--&---&-|
        # lcr w regionie funkcyjnym
        if region.begin != region.end:
            if lcr.begin >= region.begin and lcr.end <= region.end:
                region_cover = 100 * (lcr.end - lcr.begin) / (region.end - region.begin)
                return 1, region_cover

                # ----|--&-----|-&----
            # lcr zahacza początkiem o region
            elif lcr.begin > region.begin and lcr.begin < region.end:
                lcr_cover = 100 * (region.end - lcr.begin) / (lcr.end - lcr.begin)
                region_cover = 100 * (region.end - lcr.begin) / (region.end - region.begin)
                return lcr_cover, region_cover

            # --&--|--&-----|-----
            elif region.begin < lcr.end and region.end > lcr.begin:
                lcr_cover = 100 * (lcr.end - region.begin) / (lcr.end - lcr.begin)
                region_cover = 100 * (lcr.end - region.begin) / (region.end - region.begin)
                return lcr_cover, region_cover

            # --&--|-------|---&--
            elif region.begin < lcr.begin and region.end > lcr.end:
                region_cover = 1
                lcr_cover = 100 * (region.end - region.begin) / (lcr.end - lcr.begin)
                return lcr_cover, region_cover
            else:
                return 0, 0
        else:
            return 0, 0


if __name__ == "__main__":
    ldata = LoadDatabase("./data/input/SARS-2_vs_human_lcrblast_output_csv.csv", "./data/databases/uniprot/")
    ldata.read_input()

    elm_database = ELMDatabase(ldata, "./data/databases/elm/elm_instances.tsv")
    elm_database.prepare_elm_data("./data/databases/elm/elm_instances.tsv")
    ldata.add_region_info(elm_database.parsed_data)

    uniprot_db = UniprotData(ldata, "./data/databases/uniprot/")
    uniprot_db.parse_files()
    # uniprot_db.prepare_uniprot_data_to_download()
    ldata.add_region_info(uniprot_db.parsed_data)

    prosite_db = PrositeData(ldata, "./data/databases/prosite/")
    prosite_db.parse_files()
    # prosite_db.prepare_prosite_data_to_download()
    ldata.add_region_info(prosite_db.parsed_data)

    cdd_db = CDDData(ldata, "./data/databases/cd-search/")
    # cdd_db.prepare_cdd_data_to_download()
    cdd_db.parse_files()
    ldata.add_region_info(cdd_db.parsed_data)

    interpro = InterproData(ldata, "./data/databases/interpro/")
    # interpro.prepare_interpro_data_to_download()
    interpro.parse_files()
    ldata.add_region_info(interpro.parsed_data)

    print(ldata.regions_info)
    ldata.check_covarage_with_lcr()
    ldata.make_summary("./stats_SARS-2_vs_human_lcrblast_output_csv.csv")

    ldata = LoadDatabase("./data/input/SARS-2 vs human motiflcr_lcrblast_output_e001_csv.csv",
                         "./data/databases/uniprot/")
    ldata.read_input()

    elm_database = ELMDatabase(ldata, "./data/databases/elm/elm_instances.tsv")
    elm_database.prepare_elm_data("./data/databases/elm/elm_instances.tsv")
    ldata.add_region_info(elm_database.parsed_data)

    uniprot_db = UniprotData(ldata, "./data/databases/uniprot/")
    # uniprot_db.prepare_uniprot_data_to_download()
    uniprot_db.parse_files()
    ldata.add_region_info(uniprot_db.parsed_data)

    prosite_db = PrositeData(ldata, "./data/databases/prosite/")
    # prosite_db.prepare_prosite_data_to_download()
    prosite_db.parse_files()
    ldata.add_region_info(prosite_db.parsed_data)

    cdd_db = CDDData(ldata, "./data/databases/cd-search/")
    # cdd_db.prepare_cdd_data_to_download()
    cdd_db.parse_files()
    ldata.add_region_info(cdd_db.parsed_data)

    interpro = InterproData(ldata, "./data/databases/interpro/")
    # interpro.prepare_interpro_data_to_download()
    interpro.parse_files()
    ldata.add_region_info(interpro.parsed_data)

    print(ldata.regions_info)
    ldata.check_covarage_with_lcr()
    ldata.make_summary("./stats_SARS-2 vs human motiflcr_lcrblast_output_e001_csv.csv")

    ldata = LoadDatabase("./data/input/gbsc_results_nsps_2.csv",
                         "./data/databases/uniprot/")
    ldata.read_gbsc_input()

    elm_database = ELMDatabase(ldata, "./data/databases/elm/elm_instances.tsv")
    elm_database.prepare_elm_data("./data/databases/elm/elm_instances.tsv")
    ldata.add_region_info(elm_database.parsed_data)
    #
    uniprot_db = UniprotData(ldata, "./data/databases/uniprot/")
    # uniprot_db.prepare_uniprot_data_to_download()
    uniprot_db.parse_files()
    ldata.add_region_info(uniprot_db.parsed_data)

    prosite_db = PrositeData(ldata, "./data/databases/prosite/")
    # prosite_db.prepare_prosite_data_to_download()
    prosite_db.parse_files()
    ldata.add_region_info(prosite_db.parsed_data)

    cdd_db = CDDData(ldata, "./data/databases/cd-search/")
    # cdd_db.prepare_cdd_data_to_download()
    cdd_db.parse_files()
    ldata.add_region_info(cdd_db.parsed_data)

    interpro = InterproData(ldata, "./data/databases/interpro/")
    # interpro.prepare_interpro_data_to_download()
    interpro.parse_files()
    ldata.add_region_info(interpro.parsed_data)

    print(ldata.regions_info)
    ldata.check_covarage_with_lcr()
    ldata.make_summary("./stats_gbsc_results_nsps_2.csv")

    ldata = LoadDatabase("./data/input/gbsc_results_nsps.csv",
                         "./data/databases/uniprot/")
    ldata.read_gbsc_input()
    #
    elm_database = ELMDatabase(ldata, "./data/databases/elm/elm_instances.tsv")
    elm_database.prepare_elm_data("./data/databases/elm/elm_instances.tsv")
    ldata.add_region_info(elm_database.parsed_data)

    uniprot_db = UniprotData(ldata, "./data/databases/uniprot/")
    # uniprot_db.prepare_uniprot_data_to_download()
    uniprot_db.parse_files()
    ldata.add_region_info(uniprot_db.parsed_data)

    prosite_db = PrositeData(ldata, "./data/databases/prosite/")
    # prosite_db.prepare_prosite_data_to_download()
    prosite_db.parse_files()
    ldata.add_region_info(prosite_db.parsed_data)

    cdd_db = CDDData(ldata, "./data/databases/cd-search/")
    # cdd_db.prepare_cdd_data_to_download()
    cdd_db.parse_files()
    ldata.add_region_info(cdd_db.parsed_data)

    interpro = InterproData(ldata, "./data/databases/interpro/")
    # interpro.prepare_interpro_data_to_download()
    interpro.parse_files()
    ldata.add_region_info(interpro.parsed_data)

    print(ldata.regions_info)
    ldata.check_covarage_with_lcr()
    ldata.make_summary("./stats_gbsc_results_nsps.csv")
