## Variability & Dimensionality Reduction Analyst

**Last Updated:** November 29, 2024
**Your Tasks:** Run correlation analyses for Figures 3a, 3b, and 3f, then assemble the final Figure 3

---

## **WHAT Variability & Dimensionality Reduction Analyst ROLE IS ABOUT**

analyzing how **gene expression variability** changes during **yeast aging**. Specifically:

- **Figure 3a:** Do 16-hour old cells with more cell divisions express more genes?
- **Figure 3b:** Do 36-hour old cells with more cell divisions express more genes?
- **Figure 3f:** Does FIT3 gene expression decrease as cells divide more?

**The Good News:** Anaita already wrote the R scripts! just need to run them and understand the results.

---

## **STEP 1: Open RStudio (Easier for Beginners)**

Instead of running R in the terminal, let's use RStudio which is more visual:

1. Open **RStudio** application on your Mac
2. Go to **File → Open Project**
3. Navigate to `/Users/brodeur.k/Desktop/BINF6310-Group-Project`
4. Open `BINF6310-Group-Project.Rproj`

This automatically sets your working directory correctly!

---

## **STEP 2: Understand What Scripts You Have**

In RStudio, look at the **Files** pane (bottom-right). Navigate to the `scripts/` folder:

### **Scripts Anaita Created:**

1. **`figure3a_16h_correlation.R`**
   - Analyzes 16-hour old cells
   - Creates a scatter plot: Generation vs Genes Detected
   - Calculates correlation (R value and p-value)

2. **`figure3b_36h_correlation.R`**
   - Analyzes 36-hour old cells
   - Same type of plot as 3a, but for older cells

3. **Need to find:** `figure3f_fit3_correlation.R` (Anaita mentioned it exists)
   - Analyzes FIT3 gene specifically
   - Shows how FIT3 expression changes with cell divisions

---

## **STEP 3: Run Figure 3a Script**

### **In RStudio:**

1. Click on `scripts/figure3a_16h_correlation.R` to open it
2. Read through the comments (lines starting with #) to understand what it does
3. Click **Source** button (top-right of script pane) OR press **Cmd+Shift+S**

### **What Should Happen:**

You'll see output in the **Console** (bottom-left) like:
```
Total cells loaded: 125
Cells per subgroup:
16h/F 16h/S  2h  36h/F 36h/S
   ##    ##   ##    ##     ##

16h cells: XX

16h correlation: R = 0.XX, p = X.XXe-XX
Paper reported: R = 0.62, p = 1.6e-05

✓ Figure 3a saved successfully!
  - figures/figure3a_correlation_16h.pdf
  - figures/figure3a_correlation_16h.png
```

### **Check Your Figure:**

In the **Files** pane, navigate to `figures/` folder and click on:
- `figure3a_correlation_16h.png` to view your plot!

**Expected:** A scatter plot with red dots and a blue trend line.

---

## **STEP 4: Run Figure 3b Script**

Same process as Step 3, but with `figure3b_36h_correlation.R`:

1. Open the script
2. Click **Source**
3. Check the output
4. View your figure in `figures/figure3b_correlation_36h.png`

---

## **STEP 5: Find and Run Figure 3f Script**

Let's search for the FIT3 script. In RStudio **Console**, type:

```r
# List all files containing "fit" or "3f"
list.files("scripts/", pattern = "fit|3f", ignore.case = TRUE)
```

**If it exists:** Run it the same way as 3a and 3b.

**If it doesn't exist:** We'll need to create it (I can help you with this!)

---

## **STEP 6: Compare Your Results to the Paper**

Create a simple comparison table. In RStudio, create a new script called `kyle_results_summary.R`:

```r
# Kyle's Results Summary
# Comparing my results to the paper

cat("=== CORRELATION RESULTS COMPARISON ===\n\n")

cat("Figure 3a (16h):\n")
cat("  My result:    R = [your R value], p = [your p value]\n")
cat("  Paper result: R = 0.62, p = 1.6e-05\n")
cat("  Match? [YES/NO - if R values are within 0.1]\n\n")

cat("Figure 3b (36h):\n")
cat("  My result:    R = [your R value], p = [your p value]\n")
cat("  Paper result: R = 0.57, p = 2.6e-06\n")
cat("  Match? [YES/NO]\n\n")

cat("Figure 3f (FIT3 at 16h):\n")
cat("  My result:    R = [your R value], p = [your p value]\n")
cat("  Paper result: R = -0.55, p = 1.3e-04\n")
cat("  Match? [YES/NO]\n\n")

cat("Figure 3f (FIT3 at 36h):\n")
cat("  My result:    R = [your R value], p = [your p value]\n")
cat("  Paper result: R = -0.62, p = 5.6e-06\n")
cat("  Match? [YES/NO]\n\n")
```

Fill in the `[your values]` with what you got from running the scripts.

---

## **STEP 7: Assemble Final Figure 3**

Once all individual figures are created, you'll combine them into one big figure.

---

## **NEED HELP?**

### **If a script gives an error about a missing package:**

In RStudio Console, run:
```r
install.packages("package_name")
# For example:
install.packages("ggplot2")
install.packages("dplyr")
```

### **If values don't match the paper:**

This is OKAY! From the group chat, Anaita found:
- FIT3 correlations matched well (✓)
- Genes-detected correlations were slightly different
- Reason: Different DESeq2 version, slightly different normalization

explain this in your presentation!

---

## **updates**

After running these scripts,

Results:
- Figure 3a (16h): R = [your value] (paper: 0.62)
- Figure 3b (36h): R = [your value] (paper: 0.57)
- Figure 3f FIT3: R = [your values] (paper: -0.55 and -0.62)

[Say if they match or not]


## **YOUR NEXT STEPS (After reading this guide):**

1. [ ] Open RStudio and open the project
2. [ ] Run figure3a script
3. [ ] Run figure3b script
4. [ ] Find/run figure3f script
5. [ ] Fill in your results summary
6. [ ] Post update to group chat
7. [ ] create the final assembled figure

---
