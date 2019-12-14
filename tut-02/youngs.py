# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 21:15:00 2018

@author: HuangMing
"""
import argparse
import matplotlib
matplotlib.use('Agg')
import sys
import matplotlib.pyplot as plt
import numpy as np
import os.path

def rmLine(lst, key): #remove the line start with ';'
    new_lst = lst.copy()
    for item in lst:
        if item.startswith(key, 0, 5):
            new_lst.remove(item)
    return new_lst

def ReadFile(fileName):
    fileName = fileName + '.xvg'
    if not os.path.isfile(fileName):
        return False
    
    f = open(fileName, 'r')
    line = f.read().split('\n')
    line = rmLine(line, '#')
    line = rmLine(line, '@')
    f.close()
    l0 = float(line[0].split()[1]) # change according to the selection
    f0 = float(line[0].split()[2])
    strain = [] 
    stress = []
    for l in line:
        cont = l.split()
        if cont != []:
            strain.append((float(cont[1])-l0)/l0)
            stress.append(-float(cont[2])*0.0001)
    return strain, stress

def MergeData(fileName, num):
    strain_tot = []
    stress_tot = []
    for i in range(num):
        name = fileName.replace('XXX', str(i))
        print(name)
        if not ReadFile(name):
            print(name, 'not find')
            continue
        strain, stress = ReadFile(name)
        for st in strain:
            strain_tot.append(st)
        for ss in stress:
            stress_tot.append(ss)
    return strain_tot, stress_tot

def AnaData(strain, stress):
    avg_l_tmp = []
    avg_p_tmp = []
    
    avg_l = []
    avg_p = []
    p_tmp = []
    strain_tot = []
    stress_tot = []
    i = 0
    step = np.arange(0, 0.2, 0.002)
    for i in range(len(step)):
        avg_l.append([])
        avg_p.append([])
        p_tmp.append([])
    for i in range(len(step)):
        if i > 0:
            step1 = [step[i-1], step[i]]
            for ii in range(len(strain)):
                if float(step1[0]) < float(strain[ii]) and float(strain[ii]) < float(step1[1]):
                    avg_l[i-1].append(float(strain[ii]))
                    avg_p[i-1].append(float(stress[ii]))
    strain_tot = avg_l.copy()
    stress_tot = avg_p.copy()
    
    for i in range(len(avg_l)):
        if len(avg_p[i]) > 0:
            avg_l[i] = np.average(avg_l[i])
            avg_p[i] = [np.std(avg_p[i])/np.sqrt(len(avg_p)), np.average(avg_p[i])]
            p_tmp[i] = np.average(avg_p[i])
        else:
            avg_p[i] = [[], []]
            continue
    # remove empty cell
    avg_l_tmp = avg_l.copy()
    avg_p_tmp = avg_p.copy()
    
    avg_l = []
    avg_p = []
    pp = []
    for i in range(len(avg_l_tmp)):
        if avg_l_tmp[i] == []:
            continue
        else:
            avg_l.append(avg_l_tmp[i])
            avg_p.append(avg_p_tmp[i])
            pp.append(p_tmp[i])
    # Find lowest point
    minStress = min(pp[:10])
    idx = pp.index(minStress)
    
    strain_crt = avg_l[idx]
    stress_crt = avg_p[idx][1]
    print('strain: ', idx, 'stress: ', avg_p[idx][1])
    for i in range(len(avg_l)):
        if avg_l[i] == []:
            continue
        avg_l[i] = avg_l[i] - strain_crt
        avg_p[i] = [avg_p[i][0], avg_p[i][1] - stress_crt]
    return avg_l, avg_p, strain_tot, stress_tot

def GetYield(l, p):
    l_tmp = []
    p_tmp = []
    for i in range(len(l)):
        if l[i] < 0.2:
            l_tmp.append(l[i])
            p_tmp.append(p[i])
    yield_stress = max(p_tmp)
    idx = p_tmp.index(yield_stress)
    yield_strain = l_tmp[idx]
    return yield_strain, yield_stress

def StrainPlot(name, i):
    # Strain rate
    name = name + '.xvg'
    name = name.replace('XXX', str(i))
    print(name)
    if not os.path.isfile(name):
        return False
    print(name)
    f = open(name, 'r')
    line = f.read().split('\n')
    line = rmLine(line, '#')
    line = rmLine(line, '@')
    f.close()
    l0 = float(line[0].split()[1]) # change according to the selection

    strain = [] 
    time = []
    for l in line:
        cont = l.split()
        if cont != []:
            strain.append((float(cont[1])-l0)/l0)
            time.append(float(cont[0]))
    plt.figure()
    plt.title(name)
    
        
    fit2 = np.polyfit(time, strain, 1)
    fit2_fn = np.poly1d(fit2)
    print('strain rate: ', fit2[0])
    plt.xlabel('ps')
    plt.ylabel('Strain')
    plt.plot(time, strain, '.', time, fit2_fn(time), '--')
    plt.text(0.35, 0.12, 'y = {}x'.format(round(fit2[0],5)))
    plt.savefig('strain-{}.png'.format(name))

def errorCal(strain, stress, srate): # Use all data before average calculate the slope error
    l, p = flatLst(strain, stress) 
    # l and p contains all raw data from the simulations
    # if want to see the raw data, draw the plot in this function
    x = []
    y = []
    for idx in range(len(l[1:])):
        i = l[idx + 1]
        j = p[idx + 1]
        if i <= srate and i >= 0 and j > -0.08:
            x.append(i)
            y.append(j)
    fit, error = np.polyfit(x[:], y[:], 1, cov=True)
    slope_err = np.sqrt(error[0][0])

    return slope_err
    
def StressPlot(l, p, outName, slope_err, srate):
    l_tmp = l.copy()
    p_tmp = p.copy()
# =============================================================================
#   Extract min, max, avg from p
    avg_p = []
    std_p = []
    for i in range(len(p_tmp)):
        
        if p_tmp[i][0] == []:
            l.remove(l_tmp[i])
            p.remove(p_tmp[i])
        elif l_tmp[i] < 0:
            l.remove(l_tmp[i])
            p.remove(p_tmp[i])
        else:
            avg_p.append(p_tmp[i][1])
            std_p.append(p_tmp[i][0])
# =============================================================================
    
    
    x = []
    y = []
    std = []
    for idx in range(len(l)):
        i = l[idx]
        j = avg_p[idx]
        k = std_p[idx]
        if i <= srate and i >= 0 and j > -0.08:
            x.append(i)
            y.append(j)
            std.append(k)
    fit, error = np.polyfit(x[:], y[:], 1, cov=True)
    slope, err, fit_x, fit_y = slopeErrCal([x[1:], y[1:], std[1:]])
    print(fit)
    print('Slope: {}\tStd: {}'.format(slope, slope_err))
    plt.figure()
    plt.xlabel('Strain')
    plt.ylabel('Stress(GPA)')
    plt.ylim([0, 0.2])
    
    xlen = np.linspace(0, 0.05, 10)
    fit_fn = np.poly1d([slope, 0])
    
    # Calculate error
    plt.errorbar(l, avg_p, yerr=std_p)
    slopeMax = []
    slopeMin = []
    
    for idx in range(5):
        if idx == 0:
            continue
        idx = -idx
        slopeMax.append(y[idx]+std[idx]/x[idx])
        slopeMin.append(y[idx]-std[idx]/x[idx])
    
    
    slope_min = min(slopeMin)
    slope_max = max(slopeMax)
    
    delta_slope = (slope_max - slope_min)/2
    print('Slope error: ', delta_slope)
    
    
    plt.plot(l, avg_p, '.', xlen, fit_fn(xlen), '--', label='fit line')
    strain, stress = GetYield(l, avg_p)

    plt.text(0.05, 0.05, 'y = {}x, youngs modulus: {} +/- {}'.format(round(slope, 2), round(slope, 2), round(delta_slope, 2)))
    plt.tight_layout()
    plt.savefig(outName)
    plt.close('all')
    return round(slope, 2), round(delta_slope, 2), strain, stress

def summation(inList, n):
    if n > len(inList):
        print('Repeat step larger than the list length, please check!')
        sys.exit()
        
    outSum = 0
    for i in range(len(inList)):
        if i <= n - 1:
            outSum += inList[i]
    return outSum

def slopeErrCal(inList): # a stands for the slope
    tmp_s = []
    tmp_x = []
    tmp_y = []
    n = len(inList[0])
    for i in range(len(inList[0])):
        tmp_x.append(inList[0][i]**2/inList[2][i]**2)
        tmp_y.append(inList[0][i]*inList[1][i]/inList[2][i]**2)
#    
    slope = summation(tmp_y, n)/summation(tmp_x, n)
    for i in range(len(inList[0])):
        tmp_s.append((inList[1][i] - slope * inList[0][i])**2)
        
    s_init = np.sqrt(summation(tmp_s, n)/(len(tmp_s) - 1))
    x_err = np.sqrt(summation(tmp_x, n))
    slope_err = s_init/x_err
    
    fit_y = []
    fit_x = np.linspace(0, inList[0][-1], 20)
    for i in fit_x:
        fit_y.append(i*slope)
        
    return slope, slope_err, fit_x, fit_y

def flatLst(strain, stress):
    l = []
    p = []
    for i in range(len(strain)):
        for info in strain[i]:
            l.append(info)
        for info in stress[i]:
            p.append(info)
    return l, p
    
def main(name, conv, sysNum, srate):
    youngs = []
    youngs_err = []
    yield_strn = []
    yield_strs = []
    StrainPlot(name[0], 1)
    for n in name:
        youngsName = '{}-{}-youngs.png'.format(n, str(conv))
        strain, stress = MergeData(n, sysNum)
        
        avg_l, avg_p, strain_tot, stress_tot = AnaData(strain, stress)  
        
        error = 0
        slope, slope_err, y_strain, y_stress = StressPlot(avg_l, avg_p, youngsName, error, srate)
        youngs.append(slope)
        youngs_err.append(slope_err)
        yield_strn.append(y_strain)
        yield_strs.append(y_stress)

    youngs_err.sort()
    youngs_avg = str(round(np.average(youngs), 2))
    youngs_std_avg = str(round(youngs_err[1], 2))

    y_strain_avg = str(round(np.average(yield_strn), 2))
    y_stress_avg = str(round(np.average(yield_strs), 2))
    
    f = open('youngs.txt', 'w')
    str0 = 'T{}\n'.format(conv)
    info1 = 'youngs: ' + youngs_avg + '\t' + 'youngs-err: ' + youngs_std_avg + '\n'
    info2 = 'yield strain: ' + y_strain_avg + '\n'
    info3 = 'yield stress: ' + y_stress_avg + '\n'
    str_sep = '----------------------------------------------------------------'
    str1 = 'Different direction\'s youngs modulus: \n'
    str2 = 'Different direction\'s strain and stress: \n'
    f.write(str0)
    f.write(str(info1))
    f.write(str(info2))
    f.write(str(info3))
    f.write(str_sep)
    f.write(str1)
   
    for i in range(len(youngs)):
        str_tmp = '    youngs: ' + str(round(float(youngs[i]), 2)) + '\t\tyoungs_err: ' + str(round(float(youngs_err[i]), 2)) + '\n'
        f.write(str_tmp)
    
    f.write(str2)
    for i in range(len(yield_strn)):
        str_tmp = '    y_strain: ' + str(round(float(yield_strn[i]), 2)) + '\t\ty_stress: ' + str(round(float(yield_strs[i]), 2)) + '\n'
        f.write(str_tmp)
        
    f.close()

    print(str(info1))
    print(str(info2))
    print(str(info3))

    return youngs_avg, youngs_std_avg, y_strain_avg, y_stress_avg


parser = argparse.ArgumentParser(description = 'Start E cal')
parser.add_argument('-c', '-conv', dest='conv', default='100', type=str)
parser.add_argument('-sr', '-srate', dest='srate', default=0.03, type=float)
args = parser.parse_args()
conv = args.conv
srate= args.srate
multi = True
sysNum = 10
#conv = '600'
#srate = 0.05
etype= 'pressure'
name_multi = ['{}-simXXX-x-{}'.format(conv, etype), '{}-simXXX-y-{}'.format(conv, etype), '{}-simXXX-z-{}'.format(conv, etype)]
#name_multi = ['simXXX-{}-y-x'.format(conv), 'simXXX-{}-y-y'.format(conv), 'simXXX-{}-y-z'.format(conv)]

name = name_multi
data = main(name, conv, sysNum, srate)
print(data, file=sys.stderr)
