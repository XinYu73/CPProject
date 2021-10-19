import numpy as np
import unittest


class ElectromagneticField():
    def __init__(self,
                 x,
                 y,
                 z,
                 Mz=1,
                 f=20,
                 sigma=0.1,
                 mu=4 * np.pi * 10**(-7)):
        # !M = 1 A m**2
        # !f = 20 Khz
        self.x = x
        self.y = y
        self.z = z
        self.r = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        self.Mz = Mz
        self.omega = 2 * np.pi * f
        self.mu = mu
        self.sigma = sigma
        self.k = np.sqrt((0 + 1j) * self.omega * self.mu * self.sigma)

    def GetEz(self):
        self.Ez = 0

    def GetEx(self):
        self.Ex = (self.Mz / (4 * np.pi)) * (0 + 1j) * self.omega * self.mu * (
            self.y * np.exp((0 + 1j) * self.k * self.r) / self.r**3) * (
                (0 + 1j) * self.k * self.r - 1)
        print(f"At ({self.x},{self.y},{self.z}) Ex is {self.Ex} ")
        return self.Ex

    def GetEy(self):
        self.Ey = -(self.Mz / (4 * np.pi)) * (
            0 + 1j) * self.omega * self.mu * (self.x * np.exp(
                (0 + 1j) * self.k * self.r) / self.r**3) * (
                    (0 + 1j) * self.k * self.r - 1)
        print(f"At ({self.x},{self.y},{self.z}) Ex is {self.Ex} ")
        return self.Ex


class TestElectromagneticField(unittest.TestCase):
    def setUp(self):
        self.EMF = ElectromagneticField(0.5, 1, 1)

    def testEx(self):
        Ex = self.EMF.GetEx()
        self.assertTrue(np.abs(Ex - 1.0) < 10)


if __name__ == '__main__':
    unittest.main()
