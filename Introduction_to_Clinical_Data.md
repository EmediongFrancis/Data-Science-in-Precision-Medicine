<img src="materials/images/introduction-to-clinical-data-cover.png"/>

# Introduction to Clinical Data

`🕒 This module should take less than 1 hour to complete.`

`✍️ This notebook is written using Python.`

<div class="alert alert-block alert-info">
<h3>⌨️ Keyboard shortcut</h3>

These common shortcut could save your time going through this notebook:
- Run the current cell: **`Enter + Shift`**.
- Add a cell above the current cell: Press **`A`**.
- Add a cell below the current cell: Press **`B`**.
- Change a code cell to markdown cell: Select the cell, and then press **`M`**.
- Delete a cell: Press **`D`** twice.

Need more help with keyboard shortcut? Press **`H`** to look it up.
</div>

---

# Why do we need clincal tests?

Clinical tests use to measure levels of chemical components in body fluids and tissues. The goal of clincal tests are to obtain information about the health of a patient to aid in diagnosis, treatment, and prevention of disease.



# What are the different types of clinical laboratory tests?
Many different tests exist to detect and measure almost any type of chemical component such as blood glucose, electrolytes, enzymes, hormones, lipids (fats), other metabolic substances, and proteins in blood or urine [1].

The more common laboratory tests are `Urine test` `Blood tests` and `Tumor markers`.


`Urine test`: Urinalysis is a laboratory examination of urine for various cells and chemicals, such as red blood cells, white blood cells, infection, or excessive protein [2].

For exmample, blood in the urine (hematuria) may be the result of a benign (noncancerous) condition, but it can also indicate an infection or other problem. High levels of protein in the urine (proteinuria) may indicate a kidney or cardiovascular problem.

`Blood test`: Blood tests are offed used to check cell counts, measure various blood chemistries and markers of inflammation, and genetics [3].

`Tumor markers`: Tumor markers are substances either released by cancer cells into the blood or urine or substances created by the body in response to cancer cells [4].


# Plot using clinical data

In this module, we are going to plot the clincal data using python from a study published in Nature Medicine in 2019, `A longitudinal big data approach for precision health` [5].

This study collected blood samples from participants every three months for up to 8 years. The study team discovered more than 67 clinically actionable health discoveries and identified multiple molecular pathways associated with metabolic, cardiovascular and oncologic pathophysiology.

They also developed prediction models for insulin resistance by using omics measurements, illustrating their potential to replace burdensome tests. [5]


---

#### References

[1] https://stanfordhealthcare.org/medical-tests/l/lab-tests/types.html, Last Accessed September 2, 2022.

[2] https://stanfordhealthcare.org/medical-tests/l/lab-tests/types/urinalysis.html, Last Accessed September 2, 2022.

[3] https://stanfordhealthcare.org/medical-tests/l/lab-tests/types/blood.html, Last Accessed September 2, 2022.

[4] https://stanfordhealthcare.org/medical-tests/l/lab-tests/types/tumors.html, Last Accessed September 2, 2022.

[5] Schüssler-Fiorenza Rose, S. M., Contrepois, K., Moneghetti, K. J., Zhou, W., Mishra, T., Mataraso, S., ... & Snyder, M. P. (2019). A longitudinal big data approach for precision health. Nature medicine, 25(5), 792-804.

---

# What do the plots tells us?

Before we try to plot Figure 3a, 3b, 3c from the publication, let's try to understand what they mean. You could also reference the figures in publication here: https://www.nature.com/articles/s41591-019-0414-6/figures/3

These figures follow the laboratory tests for different study participants over time. They show that different people develop diabetes through different pathways, or with different initial signs of glucose dysregulation.

Each laboratory test **measures a specific aspect of glucose physiology**:

1. **FPG** - Fasting Plasma Glucose  

  It measures how well the body maintains a healthy blood glucose level in a stable environment (i.e.  homeostasis).  
  Values of 126 or above suggest that an individual might have diabetes.


2. **OGTT2HR** - Oral Glucose Tolerance Test

  It measures how well a person's body processes glucose when they receive a large amount of glucose (75gm) at once.  
  Glucose is then measured in the blood. Blood glucose levels of 200 or above are an indicator of diabetes.


3. **SSPG** - Steady-State Plasma Glucose

  It measures insulin resistance, or how well the body responds to insulin.  
  Blood glucose levels of 150 or above are considered insulin resistant. It means the body doesn’t respond as well to insulin, leading to higher blood glucose levels.


4. **HbA1C** - Hemoglobin A1C, or Glycated Hemoglobin

  It is another measure of how much sugar is in the blood. When there is a lot of sugar in the blood, it glycates (sticks to) hemoglobin.
  
  The more sugar in the blood, the more it sticks.  Since the red blood cells contain hemoglobin that last around three months, it is considered a measure of average blood glucose levels over the last three months. Values of 6.5 or above are indicative of diabetes.



## Figure 3a: Initial abnormality of OGTT2HR - Oral Glucose Tolerance Test

It shows a participant (ID: ZNDMXI3) first developed an OGTT in the diabetic range. Blood glucose levels of 200 or above are an indicator of diabetes. Here, we see measurement of 208 and 263.

<img src="materials/images/fig3a.png"/>

---

## Figure 3b: Initial abnormality of FPG - Fasting Plasma Glucose
It shows a participant (ID: ZNED4XZ) first had an abnormality in fasting glucose (glucose homeostasis). Values of 126 or above suggest that an individual might have diabetes. Here, we see multiple numeric values in red are above 126.

<img src="materials/images/fig3b.png"/>

---

## Figure 3c: Initial improvement followed by progression after viral infections

It shows a participant (ID:ZOZOWIT) had developed a diabetic range fasting glucose and hemoglobin A1C after a viral infection.

This participant was able to bring these labs back to the normal range through lifestyle changes.

However, after a second viral infection, they then developed elevated blood glucose levels, including fasting glucose, OGTT and HbA1C in the diabetic range.  


<img src="materials/images/fig3c.png"/>

---

# Create plots using clinical data from 3 study participants


```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
```

## Figure 3a: Initial abnormality of OGTT2HR - Oral Glucose Tolerance Test

Let's use data from participant (ID: ZNDMXI3) to plot the measurements from the four laboratory tests.

We will label the values in diabetic range using the color red to make them stand out. The diabetic range for the tests are:
- **FPG**: any value ≥ 126
- **OGTT2HR**: any value ≥ 200
- **SSPG**: any value ≥ 150
- **HbA1C**: any value ≥ 6.5

Check the plot after you run the code, and see if you can see that this participant had developed an OGTT blood glucose levels of 200 or above.


```python
df = pd.read_csv("data/41591_2019_414_MOESM5_ESM_3a.csv")
df
```


```python
################### PLOTS #############################
####### Plot the various data points

# This instantiates an empty plot
fig, ax = plt.subplots(figsize=(10, 5))


# This plots the fasting plasma glucose (FPG) data points
fpg, = ax.plot(df["Days"], df['FPG (mg/dl)'], color="#F1C40F", marker="o")

# This sets the FPG axis to range from 50 to 300
ax.set_ylim(50, 300)

# This adds labels to the Days and FPG axes
ax.set_xlabel("Days in Study", fontsize=12)
ax.set_ylabel('FPG (mg/dl)', fontsize=12)

# This plots the oral glucose tolerance test (OGTT) data points
ogtt, = ax.plot(df["Days"], df['OGTT_2HR (mg/dl)'], color="#3498DB", marker="o", markersize=10, linestyle="")

# This plots the insulin suppression test (SSPG) data points
sspg, = ax.plot(df["Days"], df['SSPG (mg/dl)'], color="red", marker="s", markersize=10, linestyle="")

# This adds a separate y-axis to the right of the plot representing the hemoglobin (HbA1C) axis
ax2 = ax.twinx()

# This plots the the hemoglobin (HbA1C) data points
a1c, = ax2.plot(df["Days"], df['HbA1C (%)'], color="lightgray", marker="o")

# This sets the HbA1C axis to range from 5.0 to 7.0
ax2.set_ylim(5.0, 7.0)

# This adds labels to the (HbA1C) axis
ax2.set_ylabel('HbA1C(%)', fontsize=12)



########################################## LEGEND ##########################################
#####  This anchors a legend to the top of the plot representing the various measures  #####
ax.legend(handles = [fpg, ogtt, sspg, a1c],
           labels  = ['FPG', "OGTT2HR", "SSPG", "HbA1C"],
           bbox_to_anchor=(0., 1.02, 1., .102), ncol=4,
           mode="expand", borderaxespad=0., fontsize=14);




################### ANNOTATIONS #############################

# This loops through the data;
       # data points in the diabetic range (SSPG ≥ 150) are annotated in red  ###################
xs = df["Days"]
ys = df["SSPG (mg/dl)"]

# zip joins x and y coordinates into pairs
for x,y in zip(xs,ys):
    if not np.isnan(y):  # ignore nan values
        ax.text(x, y+10, "{:.0f}".format(y), color="red", ha="center", fontsize=18);



# This loops through the data and adds numerical annotations to the "OGTT_2HR" data points;
       # data points in the diabetic range (OGTT_2HR ≥ 200) are annotated in red  ##########################
xs = df["Days"]
ys = df["OGTT_2HR (mg/dl)"]

for x,y in zip(xs,ys):
    color = "black"
    textsize=14
    if not np.isnan(y):  # ignore nan values
        if y >= 200:
            color = "red"
            textsize=18
        if y%1==0:  # to format values with decimal places appropriately
            ax.text(x, y+10, "{:.0f}".format(y), color=color, ha="center", fontsize=textsize)
        else:
            ax.text(x, y+10, "{}".format(y), color=color, ha="center", fontsize=textsize);



# This loops through the data and adds numerical annotations below the "FPG" data points #####################
xs = df["Days"]
ys = df["FPG (mg/dl)"]

for x,y in zip(xs,ys):

    label = "{:.0f}".format(y)

    ax.annotate(label, # this is the text
                 (x,y), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,-20), # distance from text to points (x,y) <-- place text below points
                 ha='center') # horizontal alignment can be left, right or center



# This loops through the data and adds numerical annotations to the  "HbA1C(%)" data points  #####################
xs = df["Days"]
ys = df["HbA1C (%)"]
for x,y in zip(xs,ys):

    label = "{}".format(y)

    ax2.annotate(label, # this is the text
                 (x,y), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,5), # distance from text to points (x,y) <-- place text below points
                 ha='center') # horizontal alignment can be left, right or center
```

---

## Figure 3b: Initial abnormality of FPG - Fasting Plasma Glucose
Let's use data from participant (ID: ZNED4XZ).

We will label the values in diabetic range using the color red to make them stand out. The diabetic range for the tests are:
- **FPG**: any value ≥ 126
- **OGTT2HR**: any value ≥ 200
- **SSPG**: any value ≥ 150
- **HbA1C**: any value ≥ 6.5

Check the plot after you run the code, and see if you can reach the same conclusion that the participant had an abnormality in fasting glucose (glucose homeostasis), where the value is at 126 or above.


```python
df = pd.read_csv("data/41591_2019_414_MOESM5_ESM_3b.csv")
df
```


```python
################################################# PLOTS ########################################################
####### Plot the various data points

# This instantiates an empty plot
fig, ax = plt.subplots(figsize=(10, 5))

# This plots the fasting plasma glucose (FPG) data points
fpg, = ax.plot(df["Days"], df['FPG (mg/dl)'].interpolate(), color="#F1C40F", marker="o")

# This sets the FPG axis to range from 60 to 190
ax.set_ylim(60, 190)

# This adds labels to the Days and FPG axes
ax.set_xlabel("Days in Study", fontsize=12)
ax.set_ylabel('FPG (mg/dl)', fontsize=12)

# This plots the oral glucose tolerance test (OGTT) data points
ogtt, = ax.plot(df["Days"], df['OGTT_2HR (mg/dl)'], color="#3498DB", marker="o", markersize=10, linestyle="")

# This selects no rows and is just used as mock (empty) data in order to display "SSPG" on the legend
sspg, = ax.plot(df["Days"][:0], df["Days"][:0], color="red", marker="s", markersize=10, linestyle="")

# This adds a separate y-axis to the right of the plot representing the hemoglobin (HbA1C) axis
ax2 = ax.twinx()

# This plots the the hemoglobin (HbA1C) data points
a1c, = ax2.plot(df["Days"], df['HbA1C (%)'].ffill(), color="lightgray", marker="o")

# This sets the HbA1C axis to range from 5.0 to 8.0
ax2.set_ylim(5.0, 8.0)

# This adds labels to the (HbA1C) axis
ax2.set_ylabel('HbA1C(%)', fontsize=12)



############################################# LEGEND ##########################################################
#####  This anchors a legend to the top of the plot representing the various measures  #####
ax.legend(handles = [fpg, ogtt, sspg, a1c],
           labels  = ['FPG', "OGTT2HR", "SSPG", "HbA1C"],
           bbox_to_anchor=(0., 1.02, 1., .102), ncol=4,
           mode="expand", borderaxespad=0., fontsize=14);




########################################## ANNOTATIONS ########################################################

# Annotations for "FPG (mg/dl)" plot ###################################

###### This gets various rows/data points to be displayed for "FPG (mg/dl)" from the dataset
df_fpg = df.iloc[[0, 4, 5, 12, 14, 20]]

#######  This loops through the various data points, adding text annotations ######
xs = df_fpg["Days"]
ys = df_fpg["FPG (mg/dl)"]

# zip joins x and y coordinates into pairs
for x,y in zip(xs,ys):
    ax.text(x, y+5, "{:.0f}".format(y), color="black", ha="center", fontsize=14)

# This displays the values in the diabetic range (FPG ≥ 126) as red
x = df.loc[15, "Days"]
y = df.loc[15, "FPG (mg/dl)"]
ax.text(x-100, y+5, "{}".format(y), color="red", ha="center", fontsize=18)

x = df.loc[16, "Days"]
y = df.loc[16, "FPG (mg/dl)"]
ax.text(x+30, y+5, "{:.0f}".format(y), color="red", ha="center", fontsize=18)

x = df.loc[19, "Days"]
y = df.loc[19, "FPG (mg/dl)"]
ax.text(x+40, y-12, "{:.0f}".format(y), color="red", ha="center", fontsize=18)




# Annotations for "HbA1C(%)" plot  ###################################

###### This gets various non-diabetic rows/data points to be displayed
df_ac1 = df.iloc[[1, 3, 4, 7, 13, 14, 16, 20]]

# This loops through the "HbA1C(%)" data points and adds numerical annotations   #####################
xs = df_ac1["Days"]
ys = df_ac1["HbA1C (%)"]

for x,y in zip(xs,ys):

    label = "{}".format(y)

    ax2.annotate(label, # this is the text
                 (x,y), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,-20), # distance from text to points (x,y) <-- place text below points
                 ha='center', fontsize=14) # horizontal alignment can be left, right or center


# This displays various values in the diabetic range (HbA1C ≥ 6.5) as red
df_ac1 = df.iloc[[17, 18, 19]]
xs = df_ac1["Days"]
ys = df_ac1["HbA1C (%)"]

for x,y in zip(xs,ys):

    label = "{}".format(y)

    ax2.annotate("{}".format(y), (x, y), textcoords="offset points",
             xytext=(0,-10), color="red", ha="center", fontsize=14)


# Annotations for "OGTT_2HR"  ###################################
x = df.loc[0, "Days"]
y = df.loc[0, "OGTT_2HR (mg/dl)"]
ax.text(x, y+5, "{}".format(y), color="black", ha="center", fontsize=18)

x = df.loc[15, "Days"]
y = df.loc[15, "OGTT_2HR (mg/dl)"]
ax.text(x, y-15, "{}".format(y), color="black", ha="center", fontsize=18)


############## Arrow annotation ##############
ax.annotate("Antibiotics", xy=(220, 120), xytext=(200, 150),
           arrowprops=dict(width=1, headwidth=10, headlength=15,
                           color = "black"), horizontalalignment='center', fontsize=14);
```

---

## Figure 3c: Initial improvement followed by progression after viral infections

We will use data from participant (ID:ZOZOWIT).

We will label the values in diabetic range using the color red to make them stand out. The diabetic range for the tests are:
- **FPG**: any value ≥ 126
- **OGTT2HR**: any value ≥ 200
- **SSPG**: any value ≥ 150
- **HbA1C**: any value ≥ 6.5

Check the plot after you run the code, and see whether it depicts the following:
1. Participant had developed a diabetic range fasting glucose and hemoglobin A1C after a viral infection.
2. Participant was able to bring these labs back to the normal range through lifestyle changes.
3. After a second viral infection, they then developed elevated blood glucose levels, including fasting glucose, OGTT and HbA1C in the diabetic range.


```python
df = pd.read_csv("data/41591_2019_414_MOESM5_ESM_3c.csv")
df
```


```python
################### PLOTS #############################
####### Plot the various data points

# This instantiates an empty plot
fig,ax = plt.subplots(figsize=(10, 5))


# This plots the fasting plasma glucose (FPG) data points
fpg, = ax.plot(df["Days"], df['FPG (mg/dl)'], color="#F1C40F", marker="o")

# This sets the FPG axis to range from 50 to 250
ax.set_ylim(50, 250)

# This adds labels to the Days and FPG axes
ax.set_xlabel("Days in Study", fontsize=12)
ax.set_ylabel('FPG (mg/dl)', fontsize=12)

# This plots the oral glucose tolerance test (OGTT) data points
ogtt, = ax.plot(df["Days"], df['OGTT_2HR (mg/dl)'], color="#3498DB", marker="o", markersize=10, linestyle="")

# This plots the insulin suppression test (SSPG) data points
sspg, = ax.plot(df["Days"], df['SSPG (mg/dl)'], color="red", marker="s", markersize=10, linestyle="")

# This adds a separate y-axis to the right of the plot representing the hemoglobin (HbA1C) axis
ax2 = ax.twinx()


# This plots the the hemoglobin (HbA1C) data points
a1c, = ax2.plot(df["Days"], df['A1C (%)'], color="lightgray", marker="o")

# This sets the HbA1C axis to range from 3.0 to 7.5
ax2.set_ylim(3.0, 7.5)

# This adds labels to the (HbA1C) axis
ax2.set_ylabel('HbA1C(%)', fontsize=12)



########################################## LEGEND ##########################################
#####  This anchors a legend to the top of the plot representing the various measures  #####
ax.legend(handles = [fpg, ogtt, sspg, a1c],
           labels  = ['FPG', "OGTT2HR", "SSPG", "HbA1C"],
           bbox_to_anchor=(0., 1.02, 1., .102), ncol=4,
           mode="expand", borderaxespad=0., fontsize=14);




########################################## ANNOTATIONS ########################################################

# Annotations for "SSPG" plot  ###################################

x = df.loc[6, "Days"]
y = df.loc[6, "SSPG (mg/dl)"]
ax.text(x-150, y, "{:.0f}".format(y), color="black", ha="center", fontsize=14)

x = df.loc[63, "Days"]
y = df.loc[63, "SSPG (mg/dl)"]
ax.text(x+30, y-15, "{}".format(y), color="black", ha="center", fontsize=14)


# Annotations for "OGTT_2HR" plot  ###################################
# This displays the values in the diabetic range (OGTT2HR ≥ 200) as red

x = df.loc[50, "Days"]
y = df.loc[50, "OGTT_2HR (mg/dl)"]
ax.text(x, y+10, "{:.0f}".format(y), color="red", ha="center", fontsize=18)

x = df.loc[167, "Days"]
y = df.loc[167, "OGTT_2HR (mg/dl)"]
ax.text(x-150, y-5, "{:.0f}".format(y), color="red", ha="center", fontsize=18)



# Annotations for "FPG (mg/dl)" plot ###################################

# This loops through the data and displays various values in the diabetic range (FPG ≥ 126) as red
# In other words, the FPG values for samples 0, 6, 57, and 138 are shown. Those in the diabetic range are displayed in red
xs = df.loc[[0, 6, 57, 138], "Days"]
ys = df.loc[[0, 6, 57, 138], "FPG (mg/dl)"]


for x,y in zip(xs,ys):
    color = "black"
    textsize=14
    space=70
    if x > 0:
        color = "red"
        textsize = 18
        space=120
    ax.text(x-space, y, "{:.0f}".format(y), color=color, ha="center", fontsize=textsize)



# Annotations for "HbA1C(%)" plot  ###################################
# The A1C values for samples 7, 8, 10, 13, 39, 56, 92, and 142 are shown. The values in the diabetic range (i.e., ≥ 6.5) are displayed in red.
xs = df.loc[[7, 8, 10, 13, 39, 56, 92, 142], "Days"]
ys = df.loc[[7, 8, 10, 13, 39, 56, 92, 142], "A1C (%)"]

# This loops through the data and displays various values in the diabetic range (HbA1C ≥ 6.5) as red
for x,y in zip(xs,ys):
    color = "black"
    textsize=14
    if y >= 6.5:
        color = "red"
        textsize = 18
    ax2.text(x, y+.1, "{}".format(y), color=color, ha="center", fontsize=textsize)


################################### Arrow annotations ###################################
ax.annotate("Viral infection", xy=(300, 90), xytext=(300, 60),
           arrowprops=dict(width=1, headwidth=10, headlength=15,
                           color = "black"), horizontalalignment='center', fontsize=12);

ax.annotate("Lifestyle change", xy=(350, 155), xytext=(50, 180),
           arrowprops=dict(width=1, headwidth=10, headlength=15,
                           color = "black"), horizontalalignment='left', fontsize=12);

ax.annotate("Viral infection", xy=(980, 190), xytext=(900, 225),
           arrowprops=dict(width=1, headwidth=10, headlength=15,
                           color = "black"), horizontalalignment='center', fontsize=12);
```

---

# Contributions & acknowledgement

Thank the following team to work on this module:

- **Module Content:** Antony Ross, Sophia Miryam Schussler-Fiorenza Rose
- **Engineering:** Amit Dixit
- **UX/UI Design & Illustration:** Kexin Cha
- **Video Production:** Francesca Goncalves
- **Project Management:** Amir Bahmani, Kexin Cha

---

Copyright (c) 2022 Stanford Data Ocean (SDO)

All rights reserved.
