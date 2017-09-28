from sympy import integrate, evalf, log, Symbol, exp
from matplotlib.pyplot import bar, figure, plot, show
from math import sqrt, ceil, floor, pi
from asyncio import get_event_loop


class chi2_apportionment():
    def __init__(self):
        self.ExperimData = []
        self.NumOfElem = 0
        self.NumInInterval = []
        self.SturgessFormula = 0 #SturgessFormula
        self.Intensity = 0
        self.StandardDev = 0
        self.MathExpect = 0 #MathExpect
        self.ExpStatSeries = [] #ExpStatSeries
        self.NormStatSeries = [] #NormStatSeries
        self.NumParSelectDistrib = 1
        self.ExpRomanovskyCriterion = 0
        self.NormRomanovskyCriterion = 0

    async def main(self, b, s):
        self.ExperimData = b
        self.NumOfElem = len(self.ExperimData)
        print(s, flush=True)
        for i in range(self.NumOfElem):
            self.MathExpect += self.ExperimData[i]
        self.MathExpect *= (1 / self.NumOfElem)
        self.Intensity = 1 / self.MathExpect                                                            #Лямбда
        dispersion = 1 / self.NumOfElem * sum((self.ExperimData[i] - self.MathExpect) ** 2 for i in range(self.NumOfElem))    #Дисперсия
        self.StandardDev = sqrt(dispersion)                                                                                 #Среднее Квадратическое Отклонение
        cko_x = sqrt(self.StandardDev / self.NumOfElem)

        mu3, mu4 = [], []
        for _ in range(self.NumOfElem):
            mu3.append((self.ExperimData[_] - self.MathExpect) ** 3) or mu4.append((self.ExperimData[_] - self.MathExpect) ** 4)
        mu3s, mu4s = sum(mu3)/self.NumOfElem, sum(mu4)/self.NumOfElem
        nu3, nu4 = mu3s/(self.StandardDev ** 3), mu4s/(self.StandardDev ** 4)                                                       #Коэффициенты асимметрии и эксцесса

        coefVar = self.StandardDev / self.MathExpect                                                             #Коэффициент Вариаций

        self.SturgessFormula = ceil((3.3 * log(self.NumOfElem, 10) + 1).evalf())

        xi_dx = floor(min(self.ExperimData))
        dx = ceil((ceil(max(self.ExperimData)) - min(self.ExperimData)) / self.SturgessFormula)                                                       #ДельтаХ

        self.counting(xi_dx, dx)
        h2exps, h2norms, h2exp, h2norm = self.chi2_counting()

        self.NumParSelectDistrib = 1
        self.ExpNumDegrFree = self.SturgessFormula - 1 - self.NumParSelectDistrib
        self.ExpRomanovskyCriterion = abs(h2exps - self.ExpNumDegrFree) / sqrt(2 * self.ExpNumDegrFree)

        self.NumParSelectDistrib = 2
        normNumDegrFree = self.SturgessFormula - 1 - self.NumParSelectDistrib
        self.NormRomanovskyCriterion = abs(h2norms - normNumDegrFree) / sqrt(2 * normNumDegrFree)

        peexp, penorm = self.chi2(self.ExpNumDegrFree, h2exps), self.chi2(normNumDegrFree, h2norms)

        self.printing(dispersion, cko_x, coefVar, dx, xi_dx, mu3, mu4, mu3s, mu4s, nu3, nu4, h2exps, h2exp, h2norm, h2norms, peexp, penorm)
        self.histogram(floor(min(self.ExperimData)), dx)
        print(end='', sep='')

    def chi2(self, k, chi2r):
        chi2 = [[0.00016, 0.0006, 0.0039, 0.016, 0.016, 0.064, 0.148, 0.445, 1.07, 1.64, 2.7, 3.8, 5.4, 6.6, 7.9, 9.5, 10.8],
                [0.02, 0.04, 0.103, 0.211, 0.446, 0.713, 1.386, 2.41, 3.22, 4.6, 6, 7.8, 9.2, 10.6, 12.4, 13.8],
                [0.115, 0.185, 0.352, 0.584, 1.005, 1.424, 2.366, 3.67, 4.64, 6.3, 7.8, 9.8, 11.3, 12.8, 14.8, 16.3],
                [0.3, 0.43, 0.71, 1.06, 1.65, 2.19, 3.36, 4.9, 6, 7.8, 9.5, 11.7, 13.3, 14.9, 16.9, 18.5],
                [0.55, 0.75, 1.14, 1.61, 2.34, 3, 4.35, 6.1, 7.3, 9.2, 11.1, 13.4, 15.1, 16.8, 18.9, 20.5],
                [0.87, 1.13, 1.63, 2.2, 3.07, 3.83, 5.35, 7.2, 8.6, 10.6, 12.6, 15.0, 16.8, 18.5, 20.7, 22.5],
                [1.24, 1.56, 2.17, 2.83, 3.82, 4.61, 6.35, 8.4, 9.8, 12, 14.1, 16.6, 18.5, 20.3, 22.6, 24.3],
                [1.65, 2.03, 2.73, 3.49, 4.59, 5.53, 7.34, 9.5, 11.0, 13.4, 15.5, 18.2, 20.1, 22, 24.3, 26.1],
                [2.09, 2.53, 3.32, 4.17, 5.38, 6.39, 8.34, 10.7, 12.2, 14.7, 16.9, 19.7, 21.7, 23.6, 26.1, 27.9],
                [2.56, 3.06, 3.94, 4.86, 6.18, 7.27, 9.34, 11.8, 13.4, 16, 18.3, 21.2, 23.2, 25.2, 27.7, 29.6],
                [3.1, 3.6, 4.6, 5.6, 7.0, 8.1, 10.3, 12.9, 14.6, 17.3, 19.7, 22.6, 24.7, 26.8, 29.4, 31.3],
                [3.6, 7.2, 5.2, 6.3, 7.8, 9, 11.3, 12.9, 15.8, 18.5, 21, 24.1, 26.2, 28.3, 30.9, 32.9],
                [4.1, 4.8, 5.9, 7.0, 8.6, 9.9, 12.3, 15.117, 19.8, 22.4, 25.5, 27.7, 29.8, 32.5, 34.5],
                [4.7, 5.4, 6.6, 7.8, 9.5, 10.8, 13.3, 16.2, 18.2, 21.1, 23.7, 26.9, 29.1, 31.3, 34, 36.1],
                [5.2, 6.0, 7.3, 8.5, 10.3, 11.7, 14.3, 17.3, 19.3, 22.3, 25, 28.3, 30.6, 32.8, 35.6, 37.7],
                [5.6, 6.6, 8.0, 9.3, 11.2, 12.6, 15.3, 18.4, 20.5, 23.5, 26.3, 29.6, 30, 34.3, 35.7, 38.6, 40.8],
                [6.4, 7.3, 8.7, 10.1, 12, 13.5, 16.3, 19.5, 21.6, 24.8, 27.6, 31, 33.4, 35.7, 38.6, 40.8],
                [7, 7.9, 9.4, 10.9, 12.9, 14.4, 17.3, 20.6, 22.8, 26, 28.9, 32.3, 34.8, 37.2, 30.1, 42.3],
                [7.6, 8.6, 10.1, 11.7, 13.7, 15.4, 18.3, 21.7, 23.9, 27.2, 30.1, 33.7, 36.2, 38.6, 41.6, 43.8],
                [8.3, 9.2, 10.9, 12.4, 14.6, 16.3, 19.3, 22.8, 25, 28.4, 31.4, 35, 37.6, 40, 43, 45.3]]
        tmp = [0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
        x2t, pt = [], []
        for i in range(1, len(tmp) + 1):
            if chi2r < max(chi2[k-1]):
                if chi2[k - 1][i - 1] < chi2r < chi2[k - 1][i]:
                    pt.append(tmp[i - 1])
                    pt.append(tmp[i])
                    x2t.append(chi2[k - 1][i - 1])
                    x2t.append(chi2[k - 1][i])
            else:
                return 'Данного значения нет в таблице распределения Хи-квадрат!'
        a = (min(pt) - max(pt)) / (max(x2t) - min(x2t))
        b = min(pt) - a * max(x2t)

        Pe = chi2r * a + b

        print(pt, x2t, a, b)

        return Pe

    def counting(self, xi_dx, dx):
        x = Symbol('x')
        b2, self.NumInInterval = self.ExperimData, [0] * self.SturgessFormula
        self.ExpStatSeries, self.NormStatSeries = [0] * self.SturgessFormula, [0] * self.SturgessFormula
        print('{}|{}-{}|{}|{}|{}'.format('№', 'xi', 'xi+dx', 'Ni', 'N*Pi exp', 'N*Pi norm'))
        for u in range(1, self.SturgessFormula+1):
            for j in range(self.NumOfElem):
                if xi_dx+dx*(u-1) < b2[j] <= xi_dx+dx*u:
                    self.NumInInterval[u - 1] += 1
        for i in range(1, self.SturgessFormula + 1):
            xi = xi_dx
            xi_dx += dx
            self.ExpStatSeries[i - 1] = (integrate(exp(-x / self.MathExpect) / self.MathExpect,
                                             (x, xi, xi_dx))).evalf()
            self.NormStatSeries[i - 1] = integrate(
                (1 / (self.StandardDev * sqrt(2 * pi))) * exp(-((x - self.MathExpect) ** 2) / (2 * self.StandardDev ** 2)),
                (x, xi, xi_dx)).evalf()

            if self.NumInInterval[i - 1] != 0:
                print('{} | {} - {} | {} | {} | {}'.format(i, round(xi, 1), round(xi_dx, 1), self.NumInInterval[i - 1],
                                                           round(self.ExpStatSeries[i - 1] * self.NumOfElem, 3), round(self.NormStatSeries[i - 1] * self.NumOfElem, 3)))
            else:
                print('{} | {} - {} | {} | {} | {}'.format(i, round(xi, 1), round(xi_dx, 1), self.NumInInterval[i - 1],
                                                           round(self.ExpStatSeries[i - 1] * self.NumOfElem, 3), round(self.NormStatSeries[i - 1] * self.NumOfElem, 3)))
        print(end='', sep='')

    def chi2_counting(self):
        h2exp, h2norm = [], []
        mol = False
        for kol in range(self.SturgessFormula):
            if self.NumInInterval[kol] != 0:
                if mol:
                    self.ExpStatSeries[kol] += self.ExpStatSeries[kol - 1]
                    self.NormStatSeries[kol] += self.NormStatSeries[kol - 1]
                    mol = False
                h2exp.append(((self.NumInInterval[kol] - self.NumOfElem * self.ExpStatSeries[kol]) ** 2) / (self.NumOfElem * self.ExpStatSeries[kol]))
                h2norm.append(((self.NumInInterval[kol] - self.NumOfElem * self.NormStatSeries[kol]) ** 2) / (self.NumOfElem * self.NormStatSeries[kol]))
            elif self.NumInInterval[kol] == 0:
                mol = True
        return (sum(h2exp), sum(h2norm), h2exp, h2norm)

    def histogram(self, xi, dx):
        fig = figure()
        ax1 = fig.add_subplot(111)
        ax2 = fig.add_subplot(111)
        ax3 = fig.add_subplot(111)
        ax1.set_title(u'Гистограмма эмпирической плотности распределения случайной величины Х', size=9)
        ax1.grid(True, color='#f0f0f0')
        x1 = [xi + dx * i for i in range(self.SturgessFormula)]

        x, y = [((xi + xi + dx) / 2) + dx * i for i in range(self.SturgessFormula)], [self.NumInInterval[i] / (dx * self.NumOfElem) for i in range(self.SturgessFormula)]

        ax1.bar(x, y, dx, color='orange')
        ax1.axis([xi, xi + dx * self.SturgessFormula, min(y) - min(y) * 0.1, max(y) + max(y)*0.1])
        ax1.tick_params(axis='y', which='major', labelcolor='black')
        ax1.set_xticks(x1)

        y = [self.Intensity * exp(- self.Intensity * x[i]) for i in range(self.SturgessFormula)]
        ax2.plot(x, y, label='Экспонинциальное распределение', color='red')

        y = [(1 / (self.StandardDev * sqrt(2 * pi))) * exp(-((x[i] - self.MathExpect) ** 2) / (2 * self.StandardDev ** 2))
             for i in range(self.SturgessFormula)
            ]

        ax2.tick_params(axis='x', which='major', labelcolor='black')
        ax3.plot(x, y, label='Нормальное распределение', color='green')
        show()

    def printing(self, dispersion, cko_x, coefVar, dx, xi_dx, mu3, mu4, mu3s, mu4s, nu3, nu4, h2exps, h2exp, h2norm, h2norms, peexp, penorm):
        print('Общее число измерений -', self.NumOfElem)
        print("Математическое Ожидание -", self.MathExpect)
        print("Дисперсия -", dispersion)
        print("СКО =", self.StandardDev)
        print("Выборочное среднее значение Х по выборки с объемом N =", cko_x)
        print("Коэффициент вариаций =", coefVar)
        print('Коэффициент ассиметрии =', nu3)
        print('Коэффициент эксцесса =', nu4)
        print()
        print('Обратный коэффициент масштаба =',  self.Intensity)
        print("Дельта Х =", dx)
        print('Максимальное значение =', max(self.ExperimData), '\nМинимальное значение =', min(self.ExperimData), '\nКоличество интервалов -', self.SturgessFormula)
        print('Минимальное значение с округлением =', xi_dx)
        print('Максимальное значение с округлением =', xi_dx + dx * self.SturgessFormula)
        print()
        print('Список интегралов Хи-квадрат экспоненциальное', *h2exp, sep='\n')
        print('Хи-квадрат распределение экспоненциальное =', h2exps)
        print('Критерий согласия Романовского экспоненциальное =', self.ExpRomanovskyCriterion)
        print('P экспоненциальное', peexp)
        print()
        print('Список интегралов Хи-квадрат нормальное', *h2norm, sep='\n')
        print('Хи-квадрат распределение нормальное =', h2norms)
        print('Критерий согласия Романовского нормальное =', self.NormRomanovskyCriterion)
        print('P нормальное', penorm)
        print()
        print(*self.ExpStatSeries)
        print(*self.NormStatSeries)


chi2 = chi2_apportionment()

"""
#Output's data in excel file
def data_output(b, etc):
    wb = Workbook()
    ws = wb.add_sheet('1')
    for i in range(len(b)):
        ws.write(i, b[i])
    wb.save(etc)
"""

get_event_loop().run_until_complete(chi2.main([float(91.8), float(162.6), float(29.1), float(193.2), float(120.9),
                                               float(366.3), float(19.2), float(137.1), float(39.3), float(501.9),
                                               float(191.7), float(118.7), float(23.5), float(152.9), float(154.1),
                                               float(496.2), float(60.8), float(230.5), float(101.9), float(54.7),
                                               float(370.7), float(248.6), float(44.1), float(410.4), float(268.4),
                                               float(199.5), float(283.2), float(59.3), float(109), float(87.7),
                                               float(131.8), float(159.4), float(95.6), float(21.9), float(421.1),
                                               float(317.2), float(254.7), float(139.6), float(220.1), float(322.1),
                                               float(45.4), float(187.3), float(63.2), float(129.4), float(190.8),
                                               float(312.7), float(251.3), float(397.3), float(67.3), float(301.5)],
                                              'Студент Бубнов Вариант №6.1'))
print()
get_event_loop().run_until_complete(chi2.main([float(193.2), float(9.2), float(299.7), float(26.3), float(51.8),
                                               float(2.7), float(184.1), float(14.5), float(87.4), float(5.2),
                                               float(289.1), float(6.1), float(72.8), float(6.7), float(155.7),
                                               float(25.8), float(1.3), float(8.4), float(75.7), float(12.6),
                                               float(32.5), float(71.5), float(92.7), float(38.9), float(434.7),
                                               float(49.5), float(113.4), float(7.2), float(105.2), float(28.4),
                                               float(1.8), float(10.5), float(98.4), float(3.9), float(208.5),
                                               float(5.7), float(216.1), float(503.8), float(23.4), float(209.3),
                                               float(7.9), float(203.4), float(21.3), float(89.2), float(81.3),
                                               float(17.9), float(4.5), float(34.6), float(354.2), float(43.2)],
                                              'Студент Воронцов Вариант №1.5'))

"""
#Input's data (for visualizer by PythonTutor.ru)
n = 50
print([float(input().replace(',','.')) for _ in range(n)])
"""

"""
#Input's data from excel and .txt files 
def data_input(n, etc):
    b = []
    if etc:
        try:
            rb = open_workbook('example.xls')
        except (FileNotFoundError):
            rb = open_workbook('example.xlsx')
        rs = rb.sheet_by_index(0)
        for i in range(n): b.append(float(str(rs.row_values(i)).replace(',', '.').replace('[', '').replace(']', '')))
    else:
        f = open('example.txt', 'r')
        b = f.read().replace(',','.').split()
        for i in range(n): b[i] = float(b[i])
    return b
"""

"""
#Output's data in excel file
def data_output(b, etc):
    wb = Workbook()
    ws = wb.add_sheet('1')
    for i in range(len(b)):
        ws.write(i, b[i])
    wb.save(etc)
"""
