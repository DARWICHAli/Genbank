from pip import main
import ui as ui

class Genbank:

    @classmethod
    def get_kingdom_choice(cls, mainwindow):
                print("getting kingdom choice")
                selected_kingdoms = []
                if(mainwindow.checkBox_prokaryota.isChecked()):
                        selected_kingdoms = selected_kingdoms + ["prokaryota"]
                if(mainwindow.checkBox_archaea.isChecked()):
                        selected_kingdoms = selected_kingdoms + ["archaea"]
                if(mainwindow.checkBox_bacteria.isChecked()):
                        selected_kingdoms = selected_kingdoms + ["bacteria"]
                if(mainwindow.checkBox_eukaryota.isChecked()):
                        selected_kingdoms = selected_kingdoms + ["eukaryota"]
                if(mainwindow.inputKingdom.toPlainText() != ""):
                        selected_kingdoms = selected_kingdoms + [mainwindow.inputKingdom.toPlainText()]
                return selected_kingdoms

    @classmethod
    def get_region_choice(cls, mainwindow):
            selected_regions = []
            if(mainwindow.checkBox_rrna.isChecked()):
                    selected_regions = selected_regions + ["rrna"]
            if(mainwindow.checkBox_cds.isChecked()):
                    selected_regions = selected_regions + ["cds"]
            if(mainwindow.checkBox_trna.isChecked()):
                    selected_regions = selected_regions + ["trna"]
            if(mainwindow.checkBox_centromere.isChecked()):
                    selected_regions = selected_regions + ["centromere"]
            if(mainwindow.checkBox_telomere.isChecked()):
                    selected_regions = selected_regions + ["telomere"]
            if(mainwindow.checkBox_3utr.isChecked()):
                    selected_regions = selected_regions + ["3utr"]
            if(mainwindow.checkBox_5utr.isChecked()):
                    selected_regions = selected_regions + ["5utr"]
            if(mainwindow.checkBox_mobile_element.isChecked()):
                    selected_regions = selected_regions + ["mobile element"]
            if(mainwindow.checkBox_mobile_ncrna.isChecked()):
                    selected_regions = selected_regions + ["ncrna"]
            if(mainwindow.checkBox_mobile_intron.isChecked()):
                    selected_regions = selected_regions + ["intron"]
            if(mainwindow.inputRegion.toPlainText() != ""):
                    selected_regions = selected_regions + [mainwindow.inputRegion.toPlainText()]
            return selected_regions

    @classmethod
    def start(cls, mainwindow):
            print("start")
            selected_kingdoms = cls.get_kingdom_choice(mainwindow)
            selected_regions = cls.get_region_choice(mainwindow)
            print(selected_kingdoms)
            print(selected_regions)

    @classmethod
    def pause(cls):
            print("pause")

    @classmethod
    def reset(cls):
            print("reset")
