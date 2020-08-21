# Code for Shipley, Dani et al. (2020)

This repository contains all of the code necessary for XYZT registration and segmentation of ChP blood vessel images.

## /data
This folder contains supporting code necessary for loading tiffs, writing .sbx (binary) files, writing metadata, and various other support programs.

## /registration
This folder contains the code used in the various steps of registration. To get started, run **RegistrationMasterPipeine.m**, and edit the first block of code to point to your data.

## /vessel segmentation
This folder has the code necessary for segmenting PECAM-stained images, dilating the periluminal space, and comparing manually segmented images to automatically segmented images.
