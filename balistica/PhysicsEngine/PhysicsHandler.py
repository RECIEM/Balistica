# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2021
#
# Authors: Santiago Nunez-Corrales <snunezcr@gmail.com>
#          Jose Brenes-Andre <jbrenes54@gmail.com>
import abc


class PhysicsHandler(abc.ABC):
    g = 9.7803253359

    @abc.abstractmethod
    def compute(self):
        pass

    @abc.abstractmethod
    def save_csv(self, filename):
        pass

    @abc.abstractmethod
    def maxT(self):
        pass

    @abc.abstractmethod
    def maxH(self):
        pass

    @abc.abstractmethod
    def totalR(self):
        pass

    @abc.abstractmethod
    def maxDistance(self):
        pass

    @abc.abstractmethod
    def totalT(self):
        pass

    @abc.abstractmethod
    def finalTheta(self):
        pass

    @abc.abstractmethod
    def finalV(self):
        pass