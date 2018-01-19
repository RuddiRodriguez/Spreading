#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:51:58 2018

@author: ruddirodriguez
"""
import sys
def squared(x):
    y = x * x
    return y

if __name__ == '__main__':
    x = float(sys.argv[1])
    sys.stdout.write(str(squared(x)))