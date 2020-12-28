import math

class Mathematica_1:

    def __init__(self, quality):
        self.quality=int(quality)
    
    def return_factors(self) :
        self.a=float()
        self.b=float()
        if self.quality == 17:
            self.a=16.72
            self.b=11.285
        elif self.quality == 16:
            self.a=15.736
            self.b=10.621
        elif self.quality == 15:
            self.a=14.753
            self.b=9.957
        elif self.quality == 14:
            self.a=13.769
            self.b=9.293
        elif self.quality == 13:
            self.a=12.786
            self.b=8.63
        elif self.quality == 12:
            self.a=11.802
            self.b=7.966
        elif self.quality == 11:
            self.a=10.819
            self.b=7.302
        elif self.quality == 10:
            self.a=9.835
            self.b=6.638
        elif self.quality == 9:
            self.a=8.8518
            self.b=5.9743
        elif self.quality == 8:
            self.a=7.8682
            self.b=5.3105
        elif self.quality == 7:
            self.a=6.8847
            self.b=4.6467
        elif self.quality == 6:
            self.a=5.9012
            self.b=3.9829

    def calcul_depth_min(self, quality_fin):
        self.return_factors()
        calcul=math.exp((int(quality_fin)-float(self.b))/float(self.a))
        return calcul




