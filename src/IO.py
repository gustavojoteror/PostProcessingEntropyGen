#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 14:26:06 2019
"""
######################FILE NAME: IO.py###########################
#================================================================                   
#   @author: Pulkit Aggarwal                                     |
#     :Master Student,                                           |
#    :Process and Energy Departmemt,                             |
#    :TU Delft,                                                  |
#    :The Netherlands                                            |
#                                                                |
# email : pulkit.pagg@gmail.com                                  |
# 								                                |
# Description: Reads and Writes data file.			           |
#								                                |
#================================================================

import re
import sys
import os

#---------------------------------------------------------------------------------------------#
#
# Read User Input file 'MOC_Config.in' <Sample Input file>
#

def ReadUserInput(name):
    IN = {}
    infile = open(name,'r')
    for line in infile:
      words = re.split(r'=| |%|\n|#',line)
      if not any(words[0] in s for s in ['\n','%',' ','#']):
        words = list(filter(None,words))
        IN[words[0]] = words[1]
    return IN
