import math
import sys

class Mathematica_1:

    def __init__(self, quality):
        self.quality= quality
    
    def return_factors(self) :
        a=float()
        b=float()
        if int(self.quality) == 17:
            a=16.72
            b=11.285
        elif int(self.quality) == 16:
            a=15.736
            b=10.621
        elif int(self.quality) == 15:
            a=14.753
            b=9.957
        elif int(self.quality) == 14:
            a=13.769
            b=9.293
        elif int(self.quality) == 13:
            a=12.786
            b=8.63
        elif int(self.quality) == 12:
            a=11.802
            b=7.966
        elif int(self.quality) == 11:
            a=10.819
            b=7.302
        elif int(self.quality) == 10:
            a=9.835
            b=6.638
        elif int(self.quality) == 9:
            a=8.8518
            b=5.9743
        elif int(self.quality) == 8:
            a=7.8682
            b=5.3105
        elif int(self.quality) == 7:
            a=6.8847
            b=4.6467
        elif int(self.quality) == 6:
            a=5.9012
            b=3.9829
        return a, b

    def calcul_depth_min(self, quality_fin, factor_a, factor_b):
        calcul=math.exp((int(quality_fin)-float(factor_b))/float(factor_a))
        return round(calcul)

if __name__=="__main__":
    class_qual_init=Mathematica_1(sys.argv[1])
    factor_a, factor_b = class_qual_init.return_factors()
    class_qual_init.calcul_depth_min(sys.argv[2], factor_a, factor_b)


