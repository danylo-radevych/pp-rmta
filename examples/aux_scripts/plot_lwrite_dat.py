#!/bin/python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy.optimize import least_squares, minimize, LinearConstraint
import matplotlib.tri as mtri

#Initialization====================================================
set_range_user = True
single_color = True

def get_letter(text):
  if text[-1] == 'S' or text[-1] == 's':
    return 's'
  if text[-1] == 'P' or text[-1] == 'p':
    return 'p'
  if text[-1] == 'D' or text[-1] == 'd':
    return 'd'
  if text[-1] == 'F' or text[-1] == 'f':
    return 'f'

rydberg = 13.605693122994 # eV
bohr = 0.529177210903 # A
efermi = 9.4087 / rydberg # Ry


xLim = [0., 0.5]
#yLim = [-1., 1.]
yLim = [-3., 3.]
rc = [1.5]
rmt = 3.

fntsz = 40

# colors=['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', \
#   'brown', 'orange', 'black', 'gray', 'violet']
#
# colors2=['red', 'green', 'blue', 'orange', 'brown', 'cyan', \
#   'magenta', 'yellow', 'violet']
#
# colorsb=['black', 'gray', 'brown']


cmap1 = 'tab10'
cmap2 = 'tab20'

paw = True

file_name_betar = "betar.dat"
file_name_chir = "chir.dat"
file_name_vlocionr = "vlocionr.dat"

file_name_d2udrderf = "d2udrderf.dat"
file_name_dloglderf = "dloglderf.dat"
file_name_dos_nlmrf = "dos_nlmrf.dat"
file_name_dos_nlnrf = "dos_nlnrf.dat"
file_name_dos_nlrf = "dos_nlrf.dat"
file_name_dos_nrf = "dos_nrf.dat"
file_name_duderf = "duderf.dat"
file_name_dudrrf = "dudrrf.dat"
file_name_etall1rf = "etall1rf.dat"
file_name_loglrf = "loglrf.dat"
file_name_mll1rf = "mll1rf.dat"
file_name_rvfullrf = "rvfullrf.dat"
file_name_urf = "urf.dat"
file_name_vfullrf = "vfullrf.dat"
file_name_vlocaer = "vlocaer.dat"
file_name_vlociong = "vlociong.dat"
file_name_vlocionr3d = "vlocionr3d.dat"
file_name_vlocscf00rf = "vlocscf00rf.dat"
file_name_vlocscfg3d = "vlocscfg3d.dat"
file_name_vlocscfr3d = "vlocscfr3d.dat"
file_name_vlocscr00rf = "vlocscr00rf.dat"
file_name_vlocscrg3d = "vlocscrg3d.dat"
file_name_vlocscrr3d = "vlocscrr3d.dat"
file_name_vsemilocr = "vsemilocr.dat"
file_name_vsemilocrf = "vsemilocrf.dat"

#===================================================================


def get_table(file_name, nlines):
  nrow = nlines - 1
  f = open(file_name, 'r')
  for i in range(2):
    line = f.readline()
  tmp = line.split()
  ncol = len(tmp)
  name = np.zeros((ncol), dtype=object)
  for i in range(ncol):
    name[i] = str(tmp[i])
  table = np.zeros((nrow, ncol), dtype=float)
  for i in range(nrow):
    line = f.readline()
    tmp = line.split()
    for j in range(ncol):
      table[i, j] = float(tmp[j])
  f.close()
  return name, table

def get_nlines(file_name):
  f = open(file_name, 'r')
  nlines = int(0)
  for line in f:
    nlines = nlines + 1
  nlines = nlines - 1
  return nlines



# loading data =================================================================

# chir -------------------------------------------------------------------------
nlines = get_nlines(file_name_chir)
chir_name, chir = get_table(file_name_chir, nlines)

plot_name="chir"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(chir_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))

for i in range(num_lines):
  ax.plot(chir[:, 2 * i], chir[:, 2 * i + 1], \
    marker='o', markersize=5, \
    color=colors[i], linestyle='', label="$\\chi_{" + \
      chir_name[2 * i + 1] + "}(r_0)$", linewidth=1)

#
plt.title("$\\chi(r_0)$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$\\chi(r_0)$', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()


# betar ------------------------------------------------------------------------
nlines = get_nlines(file_name_betar)
betar_name, betar = get_table(file_name_betar, nlines)

plot_name="betar"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(betar_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))


for i in range(num_lines):
  ax.plot(betar[:, 2 * i], betar[:, 2 * i + 1], \
    marker='o', markersize=5, \
    color=colors[i], linestyle='', label="$\\beta_{" + \
      betar_name[2 * i + 1] + "}(r_0)$", linewidth=1)

#
plt.title("$\\beta(r_0)$", fontsize=fntsz)
plt.xlabel('$r_0$ (bohr)', fontsize=fntsz)
plt.ylabel('$\\beta(r_0)$', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()



# vlocionr ---------------------------------------------------------------------
nlines = get_nlines(file_name_vlocionr)
vlocionr_name, vlocionr = get_table(file_name_vlocionr, nlines)

# nlines = get_nlines(file_name_vlocaer)
# vlocaer_name, vlocaer = get_table(file_name_vlocaer, nlines)

plot_name="vlocionr"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')


num_lines = int(len(vlocionr_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))

for i in range(num_lines):
  ax.plot(vlocionr[:, 2 * i], vlocionr[:, 2 * i + 1], \
    marker='o', markersize=5, \
    color=colors[i], linestyle='', label="$V_{\\text{loc,ion}}^{" + \
      vlocionr_name[2 * i + 1] + "}(r_0)$", linewidth=1)

# for i in range(int(len(vlocaer_name) / 2)):
#   ax.plot(vlocaer[:, 2 * i], vlocaer[:, 2 * i + 1], \
#     marker='o', markersize=5, \
#     color=colorsb[i], linestyle='-', label="$V_{\\text{AE,ion}}^{" + \
#       vlocaer_name[2 * i + 1] + "}(r_0)$", linewidth=1)

#
plt.title("$V_{\\text{loc,ion}}(r_0)$", fontsize=fntsz)
plt.xlabel('$r_0$ (bohr)', fontsize=fntsz)
plt.ylabel('$V(r_0)$ (Ry)', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()


# vsemilocr ---------------------------------------------------------------------
nlines = get_nlines(file_name_vsemilocr)
vsemilocr_name, vsemilocr = get_table(file_name_vsemilocr, nlines)

plot_name="vsemilocr"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')


num_lines = int(len(vsemilocr_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))

for i in range(num_lines):
  ax.plot(vsemilocr[:, 2 * i], vsemilocr[:, 2 * i + 1], \
    marker='o', markersize=5, \
    color=colors[i], linestyle='', label="$V_{\\text{SL}}^{" + \
      vsemilocr_name[2 * i + 1] + "}(r_0)$", linewidth=1)


#
plt.title("$V_{\\text{SL}}(r_0)$", fontsize=fntsz)
plt.xlabel('$r_0$ (bohr)', fontsize=fntsz)
plt.ylabel('$V(r_0)$ (Ry)', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()




# vfullrf_and_vlocscr00rf ------------------------------------------------------
nlines = get_nlines(file_name_vfullrf)
vfullrf_name, vfullrf = get_table(file_name_vfullrf, nlines)

nlines = get_nlines(file_name_vlocscr00rf)
vlocscr00rf_name, vlocscr00rf = get_table(file_name_vlocscr00rf, nlines)

plot_name="vfullrf_and_vlocscr00rf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')


num_lines = int(len(vfullrf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))
num_lines2 = int(len(vlocscr00rf_name) / 2)
cmap = plt.get_cmap(cmap2)
colors2 = cmap(np.linspace(0, 1.0, num_lines))


for i in range(num_lines):
  ax.plot(vfullrf[:, 2 * i], vfullrf[:, 2 * i + 1], \
    marker='o', markersize=5, \
    color=colors[i], linestyle='', label="$V_{tot}^{" + \
      vfullrf_name[2 * i + 1] + "}(r)$", linewidth=1)

for i in range(num_lines2):
  ax.plot(vlocscr00rf[:, 2 * i], vlocscr00rf[:, 2 * i + 1], \
    marker='o', markersize=5, \
    color=colors2[i], linestyle='', label="$V_{\\text{loc,scr}}^{\\text{" + \
      vlocscr00rf_name[2 * i + 1] + "}}(r)$", linewidth=1)

#
plt.title("$V_{\\text{tot}}(r)~&~V_{\\text{loc,scr}}(r)$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$V(r)$ (Ry)', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()


# rvfullrf ---------------------------------------------------------------------
nlines = get_nlines(file_name_rvfullrf)
rvfullrf_name, rvfullrf = get_table(file_name_rvfullrf, nlines)

plot_name="rvfullrf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(rvfullrf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))

for i in range(num_lines):
  ax.plot(rvfullrf[:, 2 * i], rvfullrf[:, 2 * i + 1], \
    marker='o', markersize=5, \
    color=colors[i], linestyle='', label="$V_{tot}^{" + \
      rvfullrf_name[2 * i + 1] + "}(r)$", linewidth=1)

#
plt.title("$rV_{\\text{tot}}(r)$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$r V(r)$', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()



# vlocscf00rf ------------------------------------------------------------------
nlines = get_nlines(file_name_vlocscf00rf)
vlocscf00rf_name, vlocscf00rf = get_table(file_name_vlocscf00rf, nlines)

plot_name="vlocscf00rf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(vlocscf00rf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))


for i in range(num_lines):
  if (True):
    ax.plot(vlocscf00rf[:, 2 * i], vlocscf00rf[:, 2 * i + 1], \
     marker='o', markersize=5, \
     color="black", linestyle='', label="$V_{loc,scf00}^{" + \
       vlocscf00rf_name[2 * i + 1] + "}(r)$", linewidth=1)

#
plt.title("$V_{\\text{loc,H+XC,00}}(r)$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$V(r)$ (Ry)', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()


# vsemilocrf -------------------------------------------------------------------
nlines = get_nlines(file_name_vsemilocrf)
vsemilocrf_name, vsemilocrf = get_table(file_name_vsemilocrf, nlines)

plot_name="vsemilocrf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(vsemilocrf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))


for i in range(num_lines):
  if (True):
    ax.plot(vsemilocrf[:, 2 * i], vsemilocrf[:, 2 * i + 1], \
     marker='o', markersize=5, \
     color=colors[i], linestyle='', label="$V_{SL}^{" + \
       vsemilocrf_name[2 * i + 1] + "}(r)$", linewidth=1)

#
plt.title("$V_{SL}(r)$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$V(r)$ (Ry)', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()



# urf --------------------------------------------------------------------------
nlines = get_nlines(file_name_urf)
urf_name, urf = get_table(file_name_urf, nlines)

plot_name="urf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(urf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))


for i in range(num_lines):
  if (True):
    ax.plot(urf[:, 2 * i], urf[:, 2 * i + 1], \
     marker='o', markersize=5, \
     color=colors[i], linestyle='', label="$u^{" + \
       urf_name[2 * i + 1] + "}(r)$", linewidth=1)

#
plt.title("$u(r, \\varepsilon_{\\text{F}}) = r R(r, \\varepsilon_{\\text{F}})$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$u(r, \\varepsilon_{\\text{F}})$', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()




# duderf -----------------------------------------------------------------------
nlines = get_nlines(file_name_duderf)
duderf_name, duderf = get_table(file_name_duderf, nlines)

plot_name="duderf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(duderf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))

for i in range(num_lines):
  if (True):
    ax.plot(duderf[:, 2 * i], duderf[:, 2 * i + 1], \
     marker='o', markersize=5, \
     color=colors[i], linestyle='', label="$u^{" + \
       duderf_name[2 * i + 1] + "}(r)$", linewidth=1)

#
plt.title("$\\frac{d u(r, \\varepsilon)}{d \\varepsilon}$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$\\frac{d u(r, \\varepsilon)}{d \\varepsilon}$', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()


# dudrrf -----------------------------------------------------------------------
nlines = get_nlines(file_name_dudrrf)
dudrrf_name, dudrrf = get_table(file_name_dudrrf, nlines)

plot_name="dudrrf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(dudrrf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))


for i in range(num_lines):
  if (True):
    ax.plot(dudrrf[:, 2 * i], dudrrf[:, 2 * i + 1], \
     marker='o', markersize=5, \
     color=colors[i], linestyle='', label="$u^{" + \
       dudrrf_name[2 * i + 1] + "}(r)$", linewidth=1)

#
plt.title("$\\frac{d u(r, \\varepsilon)}{d r}$", fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$\\frac{d u(r, \\varepsilon)}{d r}$', fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()



# d2udrderf --------------------------------------------------------------------
nlines = get_nlines(file_name_d2udrderf)
d2udrderf_name, d2udrderf = get_table(file_name_d2udrderf, nlines)

plot_name="d2udrderf"
fig, ax = plt.subplots()
# plt.axvline(rc[0], color='brown', linestyle='--', label='$r_c$')
# plt.axvline(rmt, color='violet', linestyle='--', label='$r_{MT}$')
plt.axhline(0, color='gray', linestyle='--')

num_lines = int(len(d2udrderf_name) / 2)
cmap = plt.get_cmap(cmap1)
colors = cmap(np.linspace(0, 1.0, num_lines))


for i in range(num_lines):
  if (True):
    ax.plot(d2udrderf[:, 2 * i], d2udrderf[:, 2 * i + 1], \
     marker='o', markersize=5, \
     color=colors[i], linestyle='', label="$u^{" + \
       d2udrderf_name[2 * i + 1] + "}(r)$", linewidth=1)

#
plt.title("$\\frac{d^2 u(r, \\varepsilon)}{d \\varepsilon  d r}$", \
  fontsize=fntsz)
plt.xlabel('$r$ (bohr)', fontsize=fntsz)
plt.ylabel('$\\frac{d^2 u(r, \\varepsilon)}{d \\varepsilon d r}$', \
  fontsize=fntsz)
ax.tick_params(axis='both', which='major', labelsize=fntsz)
plt.subplots_adjust(wspace=0, hspace=0)
plt.legend(loc="upper right", fontsize='large', prop={'size': fntsz/2.})
plt.xlim([0., rmt])
fig.set_size_inches(8.5, 11)
plt.savefig(plot_name + "_portrait.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_portrait.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(11, 8.5)
plt.savefig(plot_name + "_landscape.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_landscape.svg", format='svg', bbox_inches='tight')
fig.set_size_inches(20, 12)
plt.savefig(plot_name + "_hd.pdf", format='pdf', bbox_inches='tight')
plt.savefig(plot_name + "_hd.png", format='png', bbox_inches='tight')
plt.savefig(plot_name + "_hd.svg", format='svg', bbox_inches='tight')
plt.show()


