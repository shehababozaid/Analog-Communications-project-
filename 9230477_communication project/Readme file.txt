Analog Communications Project - Super-heterodyne Receiver (FDM)
Student: Shehab Abo Zaid Abdel Aal Abo Zaid
ID:      9230477
Date:    Fall 2025

=======================================================
1. SETUP & REQUIREMENTS
=======================================================
* MATLAB Version: R2025b
* Directory: Ensure "Full_code.m" and these 5 audio files are in the SAME folder:
  1. Short_BBCArabic2.wav
  2. Short_FM9090.wav
  3. Short_QuranPalestine.wav
  4. Short_RussianVoice.wav
  5. Short_SkyNewsArabia.wav

=======================================================
2. HOW TO RUN
=======================================================
1. Open "Full_code.m" in MATLAB.
2. Press "Run" (F5).
3. WAIT: **11 Figures will appear**, and Part 4 (Clean Signal) will play AUTOMATICALLY.
4. CHECK COMMAND WINDOW: After the first sound finishes, the interactive menu will appear.

=======================================================
3. INTERACTIVE AUDIO MENU
=======================================================
Type the number and press Enter to select a mode:

* Press [4] -> CLEAN SIGNAL (Part 4)
    - Demodulated BBC Arabic station (Clear audio).

* Press [5] -> INTERFERENCE (Part 5)
    - Simulation WITHOUT RF Filter.
    - Result: Two stations mixed together (Image Frequency interference).

* Press [6] -> FREQUENCY OFFSET (Part 6)
    - Simulation with 100 Hz error.
    - Result: "Beating" or "fluttering" amplitude effect.

* Press [0] -> EXIT
    - Stops the program ( you should press 0 when finished to exit )

=======================================================
4. OPTIONAL EXPERIMENT (OFFSET)
=======================================================
To test a severe error (1000 Hz):
1. Edit in "Full_code.m": Change `offset = 100;` to `offset = 1000;`.
2. Run code and Press [6].
3. Result: High-pitched, metallic/robotic distortion (similar to a robot voice or ghost).