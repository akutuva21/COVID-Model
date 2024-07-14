#!/bin/bash

# Define default parameters
filename="results.csv"
component="APO"
mtorc_scenario=-1
ifn_scenario=-1

# Prompt user for parameters
read -p "Enter the component to filter (default: APO): " user_component
if [[ -n $user_component ]]; then
    component=$user_component
fi

read -p "Enter the MTORC Scenario (default: All): " user_mtorc_scenario
if [[ -n $user_mtorc_scenario ]]; then
    if [[ $user_mtorc_scenario =~ ^[0-2]$ ]]; then
        mtorc_scenario=$user_mtorc_scenario
    else
        echo "Invalid input for MTORC Scenario. Using default value 0."
    fi
fi


read -p "Enter the IFN Scenario (default: All): " user_ifn_scenario
if [[ -n $user_ifn_scenario ]]; then
    if [[ $user_ifn_scenario =~ ^[0-2]$ ]]; then
        ifn_scenario=$user_ifn_scenario
    else
        echo "Invalid input for IFN Scenario. Using default value 0."
    fi
fi

# Run the Python script with the selected parameters
python plotting.py --filename "$filename" --component "$component" --mtorc_scenario "$mtorc_scenario" --ifn_scenario "$ifn_scenario"
