import numpy as np
import pandas as pd

file = "/net/zootopia/disk1/chaon/WORK/GxEX/environment/interactingE_standardized_imp.csv"

dfE = pd.read_csv(file)

lst = [
"eid",
"age",
"TDI",
"stair_climbing_frequency",
"moderate_activity_days_week",
"vigorous_activity_days_week",
"walking_days_week",
"time_driving",
"time_computer_use",
"time_tv",
"walking_pace",
"phone_use_length",
"computer_games",
"sleep_duration",
"morning_getting_up",
"nap_day",
"sleeplessness",
"daytime_sleeping",
"smoking_status",
"oily_fish_intake",
"non_oily_fish_intake",
"processed_meat_intake",
"poultry_intake",
"beef_intake",
"lamb_intake",
"pork_intake",
"cheese_intake",
"salt_added",
"hot_drink_temp",
"diet_variation",
"cooked_veg_intake",
"salad_veg_intake",
"fresh_fruit_intake",
"dried_fruit_intake",
"bread_intake",
"cereal_intake",
"tea_intake",
"coffee_intake",
"water_intake",
"alcohol_frequency"]

dfE = dfE.loc[:, lst].copy()

dfE.iloc[:, 1:] = (dfE.iloc[:, 1:] - np.array(dfE.iloc[:, 1:].mean(axis=0)).reshape(1, -1)) / np.array(dfE.iloc[:, 1:].std(axis=0)).reshape(1, -1)
tmp = dfE['sleep_duration'].copy() * dfE['sleep_duration'].copy()
dfE.insert(loc=lst.index('sleep_duration') + 1, column='sleep_duration2', value=tmp)

dfE.iloc[:, 1:] = (dfE.iloc[:, 1:] - np.array(dfE.iloc[:, 1:].mean(axis=0)).reshape(1, -1)) / np.array(dfE.iloc[:, 1:].std(axis=0)).reshape(1, -1)

dfE.to_csv("/net/zootopia/disk1/chaon/WORK/GxE/simuR1/data/envi.csv", index=False)
