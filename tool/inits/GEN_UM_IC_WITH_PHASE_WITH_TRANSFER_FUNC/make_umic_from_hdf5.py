#!/usr/bin/env python3.7

#############################################################################################
# This script is used for genertating the wave function needed for GAMER with ELBDM from    #
# MUSIC, included its real and imaginary parts. Be sure the requirements below are met:     #
# (1) Input configuration file (ics_example.conf here) for MUSIC is under the same folder.  #
# (2) Output generic hdf5 file (generated by the above input file)is under the same folder. #
# (3) Adjust the phidm_mass in physical constant if needed.                                 #
#                                                                                           #
# Use handle S/D to select output UM_IC as single/double precision and assign the hdf5 file #
# when running this script, e.g. ./make_umic_from_hdf5.py S(D) hdf5_filename                #
#                                                                                           #
# The code does the follwing things:                                                        #
# (1) Scan through the input file to collect parameters                                     #
# (2) Calculate the density from over-density given by MUSIC                                #
# (3) Calculate the phase from velocity by solving the poisson quation in spcetrum way:     #
#          a*phidm_mass*(div(velocity))/h_bar = del(phase)                                  #
#     div is evaluated via Richardson extrapolation with an adjustable order                #
# (4) Convert density and phase field to real and imaginary part of wave function           #
# (5) Output as binary file "UM_IC"                                                         #
#                                                                                           #
# Unit(s): MUCIS v.s. here                                                                  #
# (1) density: back ground density ; back ground density                                    #
# (2) velocity: box length*H_0	   ; m/s                                                    #
#############################################################################################

import os, h5py, subprocess, re, sys
import numpy as np

# Calculate the gradient by Richardson extrapolation (with periodic boundary condition)
def GRAD(field, axis, h, order):
    dim  = list(field.shape)
    dim.insert(0, order)
    grad = np.zeros(tuple(dim))
    for o in range(order):
        interval = 2**(order-1-o)
        grad[o]  = (np.roll(field, -interval, axis=axis) - np.roll(field, interval, axis=axis)) / (2*interval)
    for o in range(1,order):
        grad[o:] = (4.**o*grad[o:]-grad[o-1:-1]) / (4.**o-1.)
    return grad[-1]/h
# Physical Constant
phidm_mass = 8e-23 #eV
h_bar      = 1.0545718e-34
eV_to_kg   = 1.78266192162790e-36
Mpc_to_m   = 3.0856776e22
#####################################################################################################################

if len(sys.argv)!=3:
    raise RuntimeError("Input Error: Please select output UMIC as single/double precision (S/D) and file name of hdf5!")
else:
    input_file = "ics_example.conf"
    input_hdf5 = sys.argv[2]
    if sys.argv[1]   == "S":
        print("Output UMIC as single precision.")
        flag_D = False
    elif sys.argv[1] == "D":
        print("Output UMIC as double precision.")
        flag_D = True
    else:
        raise RuntimeError("Input Error: Please select output UMIC as single/double precision (S/D)!")

    filename_hdf5 = sys.argv[2]
    if os.path.isfile("./%s"%filename_hdf5):
        print("hdf5 file name extracted successfully! Filename is %s ."%filename_hdf5)
    else:
        raise RuntimeError("File %s cannot be found!"%filename_hdf5)


    output_file = "UM_IC"
    ghost_zone  = 4

    cmd         = "sed -e '/^#/d' %s | grep levelmax | awk {'print $3'}"%input_file
    level       = subprocess.check_output(cmd, shell=True).decode('utf-8')
    level       = eval(re.sub("\n","", level))
    print("Maximum level is %d ."%level)

    cmd         = "grep boxlength %s | awk '{printf $3}' "%input_file
    box_length  = subprocess.check_output(cmd, shell=True).decode('utf-8')
    box_length  = eval(re.sub("\n","", box_length))
    print("Box length is %.4f Mpc/h."%box_length)

    cmd         = "grep H0 %s | awk '{printf $3}' "%input_file
    H0          = subprocess.check_output(cmd, shell=True).decode('utf-8')
    H0          = eval(re.sub("\n","", H0))
    print("H0 is %.4f km/s/Mpc."%H0)

    cmd         = "grep zstart %s | awk '{printf $3}' "%input_file
    z           = subprocess.check_output(cmd, shell=True).decode('utf-8')
    z           = eval(re.sub("\n","", z))
    print("z_start is %.2f ."%z)

    factor      = box_length*1000.
    N           = 2**level
    h           = 1./N
    k_factor    = 2.*np.pi/N
    box_length *= Mpc_to_m/(H0/100.) # meter
    a           = 1./(1.+z)

    with h5py.File(input_hdf5,'r') as f:
        density_hdf5 = f['level_%03d_DM_rho'%level][ghost_zone:-ghost_zone, ghost_zone:-ghost_zone, ghost_zone:-ghost_zone]
        vx_hdf5      = f['level_%03d_DM_vx'%level][ghost_zone:-ghost_zone, ghost_zone:-ghost_zone, ghost_zone:-ghost_zone]
        vy_hdf5      = f['level_%03d_DM_vy'%level][ghost_zone:-ghost_zone, ghost_zone:-ghost_zone, ghost_zone:-ghost_zone]
        vz_hdf5      = f['level_%03d_DM_vz'%level][ghost_zone:-ghost_zone, ghost_zone:-ghost_zone, ghost_zone:-ghost_zone]
        f.close()

# Bug test###################################
    #growing_factor = 5./3.
    #density_hdf5 = (growing_factor*density_hdf5+1.)
    #criteria = (density_hdf5<0.)
    #print("Percentage of over-density smaller than 0.: %.8f %%."%(100*criteria.sum()/2**(3*level)) )
    #density_hdf5[criteria] = -1.
    #vx_hdf5[criteria] = 0.
    #vy_hdf5[criteria] = 0.
    #vz_hdf5[criteria] = 0.
    #density_hdf5 = (density_hdf5+1.)**0.5
    #vx_hdf5 *= factor
    #vy_hdf5 *= factor
    #vz_hdf5 *= factor
#############################################

# Normal version#############################
    criteria = (density_hdf5<-1.)
    print("Percentage of over-density smaller than -1: %.8f %%."%(100*criteria.sum()/2**(3*level)) )
    density_hdf5[criteria] = -1.
    vx_hdf5[criteria]      = 0.
    vy_hdf5[criteria]      = 0.
    vz_hdf5[criteria]      = 0.

    density_hdf5           = (density_hdf5+1.)**0.5
    vx_hdf5               *= factor
    vy_hdf5               *= factor
    vz_hdf5               *= factor
#############################################

    # Calculate div(v)
    vx_x                   = GRAD(vx_hdf5, axis = 0, h=h ,order=3)
    vy_y                   = GRAD(vy_hdf5, axis = 1, h=h ,order=3)
    vz_z                   = GRAD(vz_hdf5, axis = 2, h=h ,order=3)
    v_div                  = vx_x + vy_y + vz_z
    v_div                 *= a*phidm_mass*eV_to_kg/h_bar

    # Do forward DFT
    v_div_k                = np.fft.rfftn(v_div)
    # Do inverse Laplacian
    kx, ky, kz             = np.arange(N), np.arange(N), np.arange(N//2+1.)
    kxx, kyy, kzz          = np.meshgrid(kx, ky, kz)
    v_div_k               /= 2.*(np.cos(k_factor*kxx)+np.cos(k_factor*kyy)+np.cos(k_factor*kzz)-3.)
    v_div_k[0,0,0]         = 0.
    # Do inverse DFT
    phi_fft                = np.fft.irfftn(v_div_k)
    # Rescale to correct unit
    phi_fft               *= box_length/N**2
    phi_fft               -= phi_fft.min()

    # Convert to wave function
    wf_real                = density_hdf5*np.cos(phi_fft)
    wf_imag                = density_hdf5*np.sin(phi_fft)

    new_data_hdf5          = np.zeros((2, N, N, N))
    new_data_hdf5[0]       = np.swapaxes(wf_real,0,2)
    new_data_hdf5[1]       = np.swapaxes(wf_imag,0,2)
    if not flag_D:
        new_data_hdf5      = new_data_hdf5.astype(np.float32)
    print("Writing wave function to binary file...")
    with open(output_file,"wb") as f:
         new_data_hdf5.tofile(f)
         f.close()
