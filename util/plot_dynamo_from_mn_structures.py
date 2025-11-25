#! /usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import plot_multiple_files as plmf

dynamo2si = 2170378.6432790086 # V/m


parser = argparse.ArgumentParser(
    description=(
        "Plot dynamo with the given timestep.\n"
        "Usage example:\n"
        "    python plot_dynamo_from_mn_structures.py --step 080000"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("-st","--step", type=str, required=True,
                        help="Step string, e.g. 080000")
args = parser.parse_args()
step = args.step
# step = '080000'


# The non-constant components (m/n != 0/0) in FFT should be multipled by 2, 
# because there are two equivalent +k and -k counterparts, 
# which is missed in JOREK (jorek2_four.f90). 
# In contrast, the constant component (m/n = 0/0) is correct.

########################################################
#                                                      #
#         read the m/n = 1/1 component                 #
#                                                      #
########################################################

file_path = 'four_results/'
file_U_radial   = file_path + 'U_radial_n001_m001_step'     + step
file_B_radial   = file_path + 'B_radial_n001_m001_step'     + step
file_U_poloidal = file_path + 'U_poloidal_n001_m001_step' + step
file_B_poloidal = file_path + 'B_poloidal_n001_m001_step' + step

column_names, data = plmf.read_ascii_file(file_U_radial, 1)
psin  = data[:,0]; rhop = np.sqrt(psin)
am_Ur = data[:,1] * 2
re_Ur = data[:,2] * 2
im_Ur = data[:,3] * 2
ph_Ur = data[:,4]

column_names, data = plmf.read_ascii_file(file_B_radial, 1)
am_Br = data[:,1] * 2
re_Br = data[:,2] * 2
im_Br = data[:,3] * 2
ph_Br = data[:,4]

column_names, data = plmf.read_ascii_file(file_U_poloidal, 1)
am_Up = data[:,1] * 2
re_Up = data[:,2] * 2
im_Up = data[:,3] * 2
ph_Up = data[:,4]

column_names, data = plmf.read_ascii_file(file_B_poloidal, 1)
am_Bp = data[:,1] * 2
re_Bp = data[:,2] * 2
im_Bp = data[:,3] * 2
ph_Bp = data[:,4]

########################################################
#                                                      #
#         read the m/n = 0/0 component                 #
#                                                      #
########################################################

file_U00_radial   = file_path + 'U_radial_n000_m000_step'     + step
file_B00_radial   = file_path + 'B_radial_n000_m000_step'     + step
file_U00_poloidal = file_path + 'U_poloidal_n000_m000_step' + step
file_B00_poloidal = file_path + 'B_poloidal_n000_m000_step' + step

column_names, data = plmf.read_ascii_file(file_U00_radial, 1)
am_Ur00 = data[:,1] * 1
re_Ur00 = data[:,2] * 1
im_Ur00 = data[:,3] * 1
ph_Ur00 = data[:,4]

column_names, data = plmf.read_ascii_file(file_B00_radial, 1)
am_Br00 = data[:,1] * 1
re_Br00 = data[:,2] * 1
im_Br00 = data[:,3] * 1
ph_Br00 = data[:,4]

column_names, data = plmf.read_ascii_file(file_U00_poloidal, 1)
am_Up00 = data[:,1] * 1
re_Up00 = data[:,2] * 1
im_Up00 = data[:,3] * 1
ph_Up00 = data[:,4]

column_names, data = plmf.read_ascii_file(file_B00_poloidal, 1)
am_Bp00 = data[:,1] * 1
re_Bp00 = data[:,2] * 1
im_Bp00 = data[:,3] * 1
ph_Bp00 = data[:,4]

########################################################
#                                                      #
#         read the dynamo profiles from postproc/      #
#                                                      #
########################################################

file_exprs_ref = 'postproc/exprs_averaged_s' + step + '.dat'

column_names, data = plmf.read_ascii_file(file_exprs_ref, 1)
rhop_ref   = data[:,0]
dynamo11_ref = data[:,6]
dynamo00_ref = data[:,12]

########################################################
#                                                      #
#         calculate the dynamo with different methods  #
#                                                      #
########################################################

dynamo11_method1 = 1/2 * (re_Ur * re_Bp + im_Ur * im_Bp - re_Up * re_Br - im_Up * im_Br)
dynamo11_method2 = 1/2 * (am_Ur * am_Bp * np.cos(ph_Ur - ph_Bp) - am_Up * am_Br * np.cos(ph_Up - ph_Br))

dynamo00_method1 = 1 * (re_Ur00 * re_Bp00 - re_Up00 * re_Br00)
dynamo00_method2 = 1 * (am_Ur00 * am_Bp00 - am_Up00 * am_Br00)


########################################################
#                                                      #
#         plot the dynamo from different methods       #
#                                                      #
########################################################

plt.plot(rhop, dynamo00_method1 * dynamo2si, 'g--', label = 'dynamo00_method1')
plt.plot(rhop, dynamo00_method2 * dynamo2si, 'rx', label = 'dynamo00_method2')
plt.plot(rhop_ref, dynamo00_ref * dynamo2si, 'k--', label = 'dynamo00_ref')

plt.plot(rhop, dynamo11_method1 * dynamo2si, 'g-', label = 'dynamo11_method1')
plt.plot(rhop, dynamo11_method2 * dynamo2si, 'r+', label = 'dynamo11_method2')
plt.plot(rhop_ref, dynamo11_ref * dynamo2si, 'k-', label = 'dynamo11_ref')

# ph_Br2 = np.angle(re_Br + 1j * im_Br)
# plt.plot(ph_Br2)
# plt.plot(ph_Br, 'o')

plt.grid()
plt.xlim(0, 1)
legend = plt.legend()
legend.set_draggable(True)

plt.xlabel(r'$\rho_p$')
plt.ylabel('dynamo (V/m)')
plt.tight_layout()
plt.show()