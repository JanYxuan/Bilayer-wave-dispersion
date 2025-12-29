## Project Title
- Elastic wave dipsersion in the bilayer pre-stressed system
- Version: Dec 2025

## Description
This repository contains home-made MATLAB code for computing wave dispersion relations
in the bilayer system. 
The code accompanies the arXiv paper: https://arxiv.org/abs/2412.13262
Github link: https://github.com/JanYxuan/Bilayer-wave-dispersion.git

## File Structure
- 'MAIN_BilayerWaveDispersion.m' : User-configurable main script
- 'SolveStressedBilayerWaveDispersion.m' : Function for adaptively determining the search range of the roots
- 'SearchRoots.m' : Function for finding the local minimum of the determinants
- 'secular_equation.m' : Function for calculating the determinant of the characteristic matrix
- 'extrema.m' & 'extrema2.m' : Local minima detection function

## Usage
Configure parameters and run the main script: 'MAIN_BilayerWaveDispersion.m'

## Citation
X Feng, GY Li, Y Jiang, O Shortt-Nguyen, SH Yun. Optical coherence
	elastography measures mechanical tension in the lens and capsule. Acta
	Biomaterialia, 199:252-261, 2025.

## Contact
Yuxuan Jiang, @ Wellman Center for Photomedicine, MGH, Boston MA
Email: jiangyx96@gmail.com