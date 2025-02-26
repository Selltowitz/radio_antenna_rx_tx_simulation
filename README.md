# Radio Antenna RX/TX Simulation
This is a MATLAB project, which was made in the course "N14 Digital Broadband Communication" at HTW Berlin by Prof. Dr.-Ing. Markus Nölle and Dipl. Ing. Lutz Molle.

## Content
This project simulates a RX and TX antenna as well as the channel. It does the following:
- simulating loops for different number of antennas (nrAntennas)
- generate random Bits for sending/receiving
- assignment of the bits to modulation symbols (QPSK, 16QAM
- adding SNR (0-30dB)
- applying different channels (AWGN, Rayleigh, Rice
- implementing different combining methods: MRC, EGC, SDC or simple sumation
- error counting
- applying an error function and calculating BER
- aborting the simulation when a specific number of errors have occured
- comparing different aspects of the above and plotting them

  ## Basic how to use
  Use the main simulation functions:
  - simulation_DifferentRiceChannels.m
  - simulation_Diversity.m
  - simulation2.m
  Those are using the functions inside the folder "subfunctions"
