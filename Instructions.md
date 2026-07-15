# Soil Color Analysis — Easy Start Guide

Welcome! This guide will help you run the Soil Color Analysis project using RStudio. You don't need any coding experience—just follow these simple steps!

## 1. How the Project is Organized

Here is a simple map of what you'll find in the project folder. You only ever need to touch the **data** and **scripts** folders!

*   **`data/`**
    *   **`raw/`**: Your original Excel databases live here (`plant_color_data.xlsx` and `soil_color_database.xlsx`). Do not delete these!
    *   **`processed/`**: The scripts will automatically save intermediate math files here. You can ignore this folder.
*   **`output/`**: **This is your results folder!** All your final graphs, plots, and CSV summaries will magically appear here organized by step.
*   **`scripts/`**
    *   **`setup/`**: Contains a one-time setup script to install necessary tools.
    *   **`pipeline/`**: Contains the 7 numbered scripts you will run in order.
*   **`workflows/`**: Technical notes for developers (you don't need to touch these).

---

## 2. Step-by-Step Instructions

### Step A: Open the Project
1. Open the **RStudio** application on your computer.
2. Go to the top menu: **Session** > **Set Working Directory** > **Choose Directory...**
3. Select the main `SoilColor-Analysis` folder.

### Step B: Install the Requirements (One Time Only)
1. In the bottom-right panel of RStudio, find the **Files** tab.
2. Click into the `scripts/` folder, then `setup/`.
3. Click on `install_packages.R` to open it in your editor.
4. Look at the top-right of the code window and click the **Source** button. 
   *(This will safely download all the math tools the project needs. It might take a minute!)*

### Step C: Run the Analysis
We have 7 steps in the analysis pipeline. You will run them in order from 1 to 7. 

For each script:
1. In the **Files** tab, go to the `scripts/pipeline/` folder.
2. Click the script (starting with `01_plant_distribution.R`) to open it.
3. Click the **Source** button at the top-right of the code window.
4. Wait for it to finish (you'll see a completion message at the bottom of the screen).
5. Open the next script (e.g. `02_plant_normalization.R`), click **Source**, and repeat until you finish all 7!

---

## 3. Where Are My Results?
Once you finish running the scripts, open the **`output/`** folder on your computer! 

Inside, you'll find beautifully organized folders for each step:
*   **Images (`.png`)**: Colorful correlation heatmaps, distribution curves, and boxplots ready to be put into a presentation or paper.
*   **Spreadsheets (`.csv`)**: Data summaries showing exact p-values, which color variables are correlated, and statistical test results.

That's it! You've successfully analyzed your soil and soil color data!
