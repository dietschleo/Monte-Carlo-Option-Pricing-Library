@echo off
cd /d "%~dp0.."
g++ app/main.cpp src/utils/functions.cpp src/core/Simulation.cpp src/core/RandomNumber.cpp src/options/EuropeanOption.cpp src/options/CompoundOption.cpp src/options/AmericanOption.cpp src/options/AsianOption.cpp -o build/build.exe -mconsole -I src
pause
