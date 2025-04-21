# reward_perception
Code for reproducing the simulations and figures from the paper “Rewarding animals based on their subjective percepts is enabled by online Bayesian estimation of perceptual biases.”

# Author summary
We tackle a fundamental problem in animal training for neuroscience: how to incentivize animals to report what they truly perceive, when those perceptions are hidden from us. Our Bayesian method estimates subjective percepts and aligns rewards with internal experiences in real time.

INSTALLATION AND USAGE INSTRUCTIONS
===================================

1. **Install MatlabStan**
   - URL: https://github.com/brian-lau/MatlabStan
   - Before installing MatlabStan, you must first install CmdStan:
     https://mc-stan.org/docs/cmdstan-guide/installation.html
   - After installing CmdStan, *build* it using the appropriate make command, as described in the CmdStan guide.
   - Confirm that CmdStan works by running an example model from your **command-line terminal** (outside of MATLAB).

2. **Install Psignifit for MATLAB**
   - You can download it from either of these locations:
     https://psignifit.sourceforge.net/WELCOME.html
     or
     https://github.com/wichmann-lab/psignifit
   - After downloading, make sure to **add the Psignifit folder to your MATLAB path**.

3. **Compile the Stan Model**
   - Before running any scripts that use Stan, you need to compile the model on your machine by running:
     stanmodel_compile.m
   - Make sure the compilation completes **without errors**.
   - If you encounter problems, try debugging using ChatGPT or any other large language model (LLM).

4. **Plotting Figures (Without Running Stan)**
   - Even if MatlabStan is not working, you can still generate the figures from the paper using the scripts in the `code_plot_figures` folder.
   - These scripts load data from the `analysis_data` and `simulation_data` folders.

5. **Generating Data for Figures (Requires Working MatlabStan)**
   - If MatlabStan is working properly, you can generate the data for the figures by running the scripts in the `code_generating_data` folder.

   Scripts for each figure:

   - Figure 3 and S3: SimulationFigure3andS3.m
   - Figure 3GH: SimulationFigure3GH.m
   - Figure 5AB: SimulationFigure5AB.m
   - Figure 5CE: SimulationFigure5CE.m
   - Figure 5DF: SimulationFigure5DF.m
   - Figure 6: DataAnalysisFigure6.m
   - Figure 7A: SimulationFigure7A.m
   - Figures 7BCD and S2: DataAnalysisFigure7BCDandS2.m
   - Figure S5: SimulationFigureS5.m
   - Figure S6A: first run SimulationFigure5AB.m, then run SimulationFigureS6A.m
   - Figure S6BC: SimulationFigureS6BC.m

