#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 00:39:35 2018

@author: ruddirodriguez
"""
import math
def ten_renormalization(sigmai,kappa,betat,r0_ini,R_ini,pos):
    if (sigmai==0):
        sigma=sigmai;
    if len(pos)>1 and sigmai!=0:
        beta_ini = ((4*3.14*kappa)*betat)*(r0_ini/( R_ini**2));           
        sigma = sigmai*math.exp(beta_ini*(pos[-1]*1e-6));
    elif len(pos)==1 and sigmai!=0:
        beta_ini = ((4*3.14*kappa)*betat)*(r0_ini/( R_ini**2));           
        sigma = sigmai*math.exp(beta_ini*(pos*1e-6));
                  
    return (sigma)
    