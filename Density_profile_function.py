# this file is written for function
import os

import numpy as np
import pandas as pd
import scipy as sp

np.set_printoptions(threshold=np.inf)


def current_path():
    """
    path fucntion to confirm correct path
    """
    return print(f"Current working directory: \n {os.getcwd()} \n")


def file_read(pwd):
    """
    pwd -> the current working path
    use os.chdir to set up working directory
    pd.read_csv will skip 23 rows and read the data.xvg since then (delimiter = '\t')

    for the xy box file, will have slightly differnt
    so set up boolean condition for it
    """
    os.chdir(str(pwd))
    current_path()
    file_name = []
    file_dict = {}
    for _, file in enumerate(os.listdir()):
        temp_file_name = str(file)[:-4]
        temp_file_df = pd.read_csv(f"{file}", skiprows=23, delimiter="\t")
        temp_file_df = temp_file_df.iloc[:, 0].str.split(expand=True)
        file_name.append(temp_file_name)

        if temp_file_name.endswith("xy"):
            temp_file_df = temp_file_df[1:]
            file_dict[temp_file_name] = temp_file_df.astype(float)
        else:
            file_dict[temp_file_name] = temp_file_df.astype(float)

    return file_name, file_dict


# set up interpolate function for data
def data_interp(
    filelist,
    filedict,
    cosolvent,
    atom_cosolvent,
):
    """
    this function is defined as the interpolate function
    filelist -> the list obtained from file_read
    filedict -> the dictionary data obtained from file_read

    each data will be interpolate to 500 data points by cubic method

    co-solvent will be dived by the atom numbers to make the density profile becomae per atom

    xy box will be used to calcuate the box area

    return the filedict as the dictionary form
    """
    for _, file in enumerate(filelist):
        temp_file_df = filedict[file]
        # boolean for file name
        if file.endswith("along_z"):
            # boolean density profile of cosolvent
            if file.endswith(f"{cosolvent}_along_z"):
                temp_interp1d = sp.interpolate.interp1d(
                    temp_file_df.iloc[:, 0],
                    temp_file_df.iloc[:, 1],
                    kind="cubic",
                )
                xi = np.linspace(
                    min(temp_file_df.iloc[:, 0]),
                    max(temp_file_df.iloc[:, 0]),
                    500,
                )
                yi = temp_interp1d(xi) / atom_cosolvent

            else:
                temp_interp1d = sp.interpolate.interp1d(
                    temp_file_df.iloc[:, 0],
                    temp_file_df.iloc[:, 1],
                    kind="cubic",
                )
                xi = np.linspace(
                    min(temp_file_df.iloc[:, 0]),
                    max(temp_file_df.iloc[:, 0]),
                    500,
                )
                yi = temp_interp1d(xi)
            temp_df = pd.DataFrame(
                {"interp1d_x": xi, "interp1d_y": yi},
                columns=["interp1d_x", "interp1d_y"],
            )
            filedict[file] = pd.concat(
                [filedict[file], temp_df],
                ignore_index=False,
                axis=1,
            )
        # interpolation for the box xy
        elif file.endswith("_xy"):
            filedict[file] = temp_file_df.rename(
                columns={
                    0: "time",
                    1: "x",
                    2: "y",
                    3: "box_area",
                },
            )
            filedict[file]["box_area"] = filedict[file]["x"] * filedict[file]["y"]
        else:
            print(f"file {file} error here")
    return filedict


def cal_integra(data):
    """
    reset the index
    integrate the whole data set

    return the y_integration value and normalization value


    """
    data = data.reset_index(drop=True)
    dx_value = data.iloc[1]["interp1d_x"] - data.iloc[0]["interp1d_x"]
    y_integration = sp.integrate.cumulative_trapezoid(
        data["interp1d_y"], data["interp1d_x"], dx=dx_value, initial=0
    )
    y_normalize = y_integration / y_integration[-1]
    return y_integration, y_normalize


def find_nearest(data, value):
    """
    this function to find the nearest point for the z boundary
    will return index and the value
    """
    data = pd.Series(data)
    pos_index = (np.abs(data - value)).argmin()
    org_index = data.index[pos_index]

    return org_index, data[org_index]


def len_integra_main(dataset, time_range, system_type, cosolvent, lipid_perleaflet):
    """
    this function is used to do the calucation and get the position
    output the table composed by z position and integrations of :
    1. 90% water position (the boundary between bulk water and lipid headgroup, unit: nm)
    2. C21C31 position (the boundary between lipid headgroup and tails, unit: nm)
    3. P atom position (unit: nm)
    4. co-solvent tail integra (integration of co-solvent molecules in the tails region, unit: #/nm^2)
    5. co-solvent head integra (integration of co-solvent moleucles in the head region, unit: #/nm^2)
    6. co-solvent bulk integra (integration of co-solvent molecuels in the bulk solvents, unit: #/nm^2)
    7. bulk len (the len of co-solvent bulk integra, unit: #/nm^2)
    8. solvent integra in head (integration of sovlent in the head region, unit: #/nm^2)
    9. average box area (unit: nm^2)
    10. bulk water density number (unit: #/nm^3)
    11. DB (unit: nm)
    12. Dh (unit: nm)
    13. Dc (unit: nm)
    14. Dpp (unit: nm)
    15. Kp (partition coefficient)
    16. Ps (co-sovlent localization)
    17. Nw (number of water molecules in the headgroup region)
    """

    # find the parallel run number
    boxarea_result = []
    z_result = {}
    temp_dict = {}
    fill_data = np.zeros((17, 3))

    # set up return data frame
    output_df = pd.DataFrame(data=fill_data)
    output_df = output_df.rename(columns=lambda x: "r" + str(x))
    output_df.index = [
        "90%Water",
        "C21C31",
        "P",
        "Cosolvent_Tail",
        "Cosolvent_Head",
        "Cosolvent_Bulk",
        "Bulk_Length",
        "Water in Lipid Head",
        "Average Box Area",
        "Bulk water",
        "DB",
        "Dc",
        "Dh",
        "Dpp",
        "Kp",
        "Ps",
        "Nw",
    ]
    filelist = list(dataset.keys())
    for _, file in enumerate(filelist):
        # set up z position for each group in negative and postive sides
        # get the parallel run number
        r_ndx = file.index("_r") + 1
        if file.endswith("_xy"):
            time_range_index = np.where(dataset[file]["time"] == time_range)[0][0]
            box_area = dataset[file][time_range_index:].describe()
            output_df.loc["Average Box Area", file[r_ndx : r_ndx + 2]] = box_area[
                "box_area"
            ]["mean"]
            boxarea_result.append([file, box_area])
        elif file.endswith("along_z"):
            along_z_ndx = file.index("_along_z")
            key = file[r_ndx:along_z_ndx]
            z_result[f"{key}_z_neg"] = []
            z_result[f"{key}_z_pos"] = []
            temp_dict[f"{key}"] = []
            # water boundary
            if file[r_ndx + 2 :].startswith("_water"):
                # average 1 nm water molecules density
                # negative position
                water_neg_1nm = dataset[file].iloc[0]["interp1d_x"] + 1
                water_neg_1nm_ndx, _ = find_nearest(
                    dataset[file].iloc[:]["interp1d_x"],
                    water_neg_1nm,
                )
                average_water_neg_1nm = np.average(
                    dataset[file].iloc[0 : (water_neg_1nm_ndx + 1)]["interp1d_y"]
                )
                water90p_neg_ndx, _ = find_nearest(
                    dataset[file].iloc[0:250]["interp1d_y"], 0.9 * average_water_neg_1nm
                )
                z_neg_position = dataset[file].iloc[water90p_neg_ndx]["interp1d_x"]
                z_result[f"{key}_z_neg"].append(z_neg_position)

                # positive position
                water_pos_1nm = dataset[file].iloc[-1]["interp1d_x"] - 1
                water_pos_1nm_ndx, _ = find_nearest(
                    dataset[file].iloc[:]["interp1d_x"], water_pos_1nm
                )
                average_water_pos_1nm = np.average(
                    dataset[file].iloc[water_pos_1nm_ndx:]["interp1d_y"]
                )
                water90p_pos_ndx, _ = find_nearest(
                    dataset[file].iloc[250:]["interp1d_y"],
                    0.9 * average_water_pos_1nm,
                )
                z_pos_position = dataset[file].iloc[water90p_pos_ndx]["interp1d_x"]
                z_result[f"{key}_z_pos"].append(z_pos_position)

                # calculate the bulk water number density
                output_df.loc["Bulk water", file[r_ndx : r_ndx + 2]] = np.average(
                    [average_water_neg_1nm, average_water_pos_1nm]
                )

                # temp data storage
                temp_dict[f"{key}"] = dataset[file]

            # other
            else:
                # negative side
                _, y_normalize_neg = cal_integra(dataset[file].iloc[0:250])
                z_neg_ndx, _ = find_nearest(y_normalize_neg, 0.5)
                z_temp_df = dataset[file].iloc[0:250].reset_index()
                z_neg_position = z_temp_df.iloc[z_neg_ndx]["interp1d_x"]
                z_result[f"{key}_z_neg"].append(z_neg_position)

                # positive side
                _, y_normalize_pos = cal_integra(dataset[file].iloc[250:])
                z_pos_ndx, _ = find_nearest(y_normalize_pos, 0.5)
                z_temp_df = dataset[file].iloc[250:].reset_index()
                z_pos_position = z_temp_df.iloc[z_pos_ndx]["interp1d_x"]
                z_result[f"{key}_z_pos"].append(z_pos_position)

                # temp data storage
                temp_dict[f"{key}"] = dataset[file]

    # integration for cosolvent
    cosolvent_data = {}
    if system_type == "y":
        for i in range(0, len(output_df.columns)):
            for _, elem in enumerate(["pos", "neg"]):
                cosolvent_data[f"r{i}_Cosolvent_Bulk_integra_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Bulk_len_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Head_integra_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Head_len_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Tail_integra_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Tail_len_{elem}"] = []
                cosolvent_data[f"r{i}_water_head_integra_{elem}"] = []

            temp_data_cosolvent = pd.DataFrame(
                {
                    "interp1d_x": temp_dict[f"r{i}_{cosolvent}"].iloc[:]["interp1d_x"],
                    "interp1d_y": temp_dict[f"r{i}_{cosolvent}"].iloc[:]["interp1d_y"],
                }
            )
            # Bulk Cosovlent
            # negative
            cosolvent_1p2nm_neg = temp_data_cosolvent.iloc[0]["interp1d_x"] + 1.2
            z_neg_ndx, z_neg_value = find_nearest(
                temp_data_cosolvent.iloc[:250]["interp1d_x"], cosolvent_1p2nm_neg
            )
            bulk_cosolvent_integration_neg_arr, _ = cal_integra(
                temp_data_cosolvent[0 : (z_neg_ndx + 1)]
            )
            bulk_cosolvent_integration_neg_result = bulk_cosolvent_integration_neg_arr[
                -1
            ]
            bulk_len_neg = np.abs(
                z_neg_value - temp_data_cosolvent.iloc[0]["interp1d_x"]
            )
            cosolvent_data[f"r{i}_Cosolvent_Bulk_integra_neg"].append(
                bulk_cosolvent_integration_neg_result
            )
            cosolvent_data[f"r{i}_Cosolvent_Bulk_len_neg"].append(bulk_len_neg)
            # positive
            cosolvent_1p2nm_pos = temp_data_cosolvent.iloc[-1]["interp1d_x"] - 1.2
            z_pos_ndx, z_pos_value = find_nearest(
                temp_data_cosolvent.iloc[250:]["interp1d_x"], cosolvent_1p2nm_pos
            )
            bulk_cosolvent_integration_pos_arr, _ = cal_integra(
                temp_data_cosolvent[z_pos_ndx:501]
            )
            bulk_cosolvent_integration_pos_result = bulk_cosolvent_integration_pos_arr[
                -1
            ]
            bulk_len_pos = np.abs(
                temp_data_cosolvent.iloc[-1]["interp1d_x"] - z_pos_value
            )
            cosolvent_data[f"r{i}_Cosolvent_Bulk_integra_pos"].append(
                bulk_cosolvent_integration_pos_result
            )
            cosolvent_data[f"r{i}_Cosolvent_Bulk_len_pos"].append(bulk_len_pos)

            # Head
            # negative
            head_cosolvent_s_neg_ndx = temp_data_cosolvent["interp1d_x"][
                temp_data_cosolvent["interp1d_x"] == z_result[f"r{i}_water_z_neg"][0]
            ].index[0]
            head_cosolvent_e_neg_ndx = temp_data_cosolvent["interp1d_x"][
                temp_data_cosolvent["interp1d_x"] == z_result[f"r{i}_C21C31_z_neg"][0]
            ].index[0]

            head_cosolvent_integration_neg_arr, _ = cal_integra(
                temp_data_cosolvent[head_cosolvent_s_neg_ndx:head_cosolvent_e_neg_ndx]
            )
            head_cosolvent_integration_neg_result = head_cosolvent_integration_neg_arr[
                -1
            ]
            head_cosolvent_len_neg = np.abs(
                z_result[f"r{i}_water_z_neg"][0] - z_result[f"r{i}_C21C31_z_neg"][0]
            )
            cosolvent_data[f"r{i}_Cosolvent_Head_integra_neg"].append(
                head_cosolvent_integration_neg_result
            )
            cosolvent_data[f"r{i}_Cosolvent_Head_len_neg"].append(
                head_cosolvent_len_neg
            )
            # positive
            head_cosolvent_s_pos_ndx = (
                temp_data_cosolvent["interp1d_x"][
                    temp_data_cosolvent["interp1d_x"]
                    == z_result[f"r{i}_water_z_pos"][0]
                ].index[0]
                + 1
            )
            head_cosolvent_e_pos_ndx = (
                temp_data_cosolvent["interp1d_x"][
                    temp_data_cosolvent["interp1d_x"]
                    == z_result[f"r{i}_C21C31_z_pos"][0]
                ].index[0]
                + 1
            )
            head_cosolvent_integration_pos_arr, _ = cal_integra(
                temp_data_cosolvent[head_cosolvent_e_pos_ndx:head_cosolvent_s_pos_ndx]
            )
            head_cosolvent_integration_pos_result = head_cosolvent_integration_pos_arr[
                -1
            ]
            head_cosolvent_len_pos = np.abs(
                z_result[f"r{i}_water_z_pos"][0] - z_result[f"r{i}_C21C31_z_pos"][0]
            )
            cosolvent_data[f"r{i}_Cosolvent_Head_integra_pos"].append(
                head_cosolvent_integration_pos_result
            )
            cosolvent_data[f"r{i}_Cosolvent_Head_len_pos"].append(
                head_cosolvent_len_pos
            )

            # Tail
            # negative
            tail_cosolvent_s_neg_ndx = temp_data_cosolvent["interp1d_x"][
                temp_data_cosolvent["interp1d_x"] == z_result[f"r{i}_C21C31_z_neg"][0]
            ].index[0]
            tail_cosolvent_integration_neg_arr, _ = cal_integra(
                temp_data_cosolvent[tail_cosolvent_s_neg_ndx:250]
            )
            tail_cosolvent_integration_neg_result = tail_cosolvent_integration_neg_arr[
                -1
            ]
            tail_cosolvent_len_neg = np.abs(z_result[f"r{i}_C21C31_z_neg"][0])
            cosolvent_data[f"r{i}_Cosolvent_Tail_integra_neg"].append(
                tail_cosolvent_integration_neg_result
            )
            cosolvent_data[f"r{i}_Cosolvent_Tail_len_neg"].append(
                tail_cosolvent_len_neg
            )
            # positive
            tail_cosolvent_s_pos_ndx = (
                temp_data_cosolvent["interp1d_x"][
                    temp_data_cosolvent["interp1d_x"]
                    == z_result[f"r{i}_C21C31_z_pos"][0]
                ].index[0]
                + 1
            )
            tail_cosolvent_integration_pos_arr, _ = cal_integra(
                temp_data_cosolvent[250:tail_cosolvent_s_pos_ndx]
            )
            tail_cosolvent_integration_pos_result = tail_cosolvent_integration_pos_arr[
                -1
            ]
            tail_cosolvent_len_pos = np.abs(z_result[f"r{i}_C21C31_z_pos"][0])
            cosolvent_data[f"r{i}_Cosolvent_Tail_integra_pos"].append(
                tail_cosolvent_integration_pos_result
            )
            cosolvent_data[f"r{i}_Cosolvent_Tail_len_pos"].append(
                tail_cosolvent_len_pos
            )

            # water in headgroup
            temp_data_water = pd.DataFrame(
                {
                    "interp1d_x": temp_dict[f"r{i}_water"]["interp1d_x"],
                    "interp1d_y": temp_dict[f"r{i}_water"]["interp1d_y"],
                }
            )
            # negative
            head_water_integration_neg_arr, _ = cal_integra(
                temp_data_water[head_cosolvent_s_neg_ndx:head_cosolvent_e_neg_ndx]
            )
            head_water_integration_neg_result = head_water_integration_neg_arr[-1]
            cosolvent_data[f"r{i}_water_head_integra_neg"].append(
                head_water_integration_neg_result
            )
            # positive
            head_water_integration_pos_arr, _ = cal_integra(
                temp_data_water[head_cosolvent_e_pos_ndx:head_cosolvent_s_pos_ndx]
            )
            head_water_integration_pos_result = head_water_integration_pos_arr[-1]
            cosolvent_data[f"r{i}_water_head_integra_pos"].append(
                head_water_integration_pos_result
            )
    elif system_type == "n":
        for i in range(0, len(output_df.columns)):
            for _, elem in enumerate(["pos", "neg"]):
                cosolvent_data[f"r{i}_Cosolvent_Bulk_integra_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Bulk_len_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Head_integra_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Head_len_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Tail_integra_{elem}"] = []
                cosolvent_data[f"r{i}_Cosolvent_Tail_len_{elem}"] = []
                cosolvent_data[f"r{i}_water_head_integra_{elem}"] = []
            temp_data_water = pd.DataFrame(
                {
                    "interp1d_x": temp_dict[f"r{i}_water"]["interp1d_x"],
                    "interp1d_y": temp_dict[f"r{i}_water"]["interp1d_y"],
                }
            )
            # negative
            head_water_s_neg_ndx = temp_data_water["interp1d_x"][
                temp_data_water["interp1d_x"] == z_result[f"r{i}_water_z_neg"][0]
            ].index[0]
            head_water_e_neg_ndx = temp_data_water["interp1d_x"][
                temp_data_water["interp1d_x"] == z_result[f"r{i}_C21C31_z_neg"][0]
            ].index[0]
            head_water_integration_neg_arr, _ = cal_integra(
                temp_data_water[head_water_s_neg_ndx:head_water_e_neg_ndx]
            )
            head_water_integration_neg_result = head_water_integration_neg_arr[-1]
            head_water_len_neg = np.abs(
                temp_data_water.iloc[head_water_s_neg_ndx]["interp1d_x"]
                - temp_data_water.iloc[head_water_e_neg_ndx]["interp1d_x"]
            )
            cosolvent_data[f"r{i}_water_head_integra_neg"].append(
                head_water_integration_neg_result
            )
            cosolvent_data[f"r{i}_Cosolvent_Head_len_neg"].append(head_water_len_neg)

            # positive
            head_water_s_pos_ndx = (
                temp_data_water["interp1d_x"][
                    temp_data_water["interp1d_x"] == z_result[f"r{i}_water_z_pos"][0]
                ].index[0]
                + 1
            )
            head_water_e_pos_ndx = (
                temp_data_water["interp1d_x"][
                    temp_data_water["interp1d_x"] == z_result[f"r{i}_C21C31_z_pos"][0]
                ].index[0]
                + 1
            )
            head_water_integration_pos_arr, _ = cal_integra(
                temp_data_water[head_water_e_pos_ndx:head_water_s_pos_ndx]
            )
            head_water_integration_pos_result = head_water_integration_pos_arr[-1]
            head_water_len_pos = np.abs(
                temp_data_water.iloc[head_water_s_pos_ndx]["interp1d_x"]
                - temp_data_water.iloc[head_water_e_neg_ndx]["interp1d_x"]
            )
            cosolvent_data[f"r{i}_water_head_integra_pos"].append(
                head_water_integration_pos_result
            )
            cosolvent_data[f"r{i}_Cosovlent_head_len_pos"].append(head_water_len_pos)

    else:
        print("unknown system_type")

    # calculate the final parameter
    for i in range(0, len(output_df.columns)):
        r_n = f"r{i}"
        output_df.loc["90%Water", r_n] = (
            np.abs(z_result[f"{r_n}_water_z_neg"][0])
            + np.abs(z_result[f"{r_n}_water_z_pos"][0])
        ) / 2
        output_df.loc["C21C31", r_n] = (
            np.abs(z_result[f"{r_n}_C21C31_z_neg"][0])
            + np.abs(z_result[f"{r_n}_C21C31_z_pos"][0])
        ) / 2
        output_df.loc["P", r_n] = (
            np.abs(z_result[f"{r_n}_P_z_neg"][0])
            + np.abs(z_result[f"{r_n}_P_z_pos"][0])
        ) / 2
        output_df.loc["Cosolvent_Tail", r_n] = (
            cosolvent_data[f"{r_n}_Cosolvent_Tail_integra_neg"][0]
            + cosolvent_data[f"{r_n}_Cosolvent_Tail_integra_pos"][0]
        ) / 2
        output_df.loc["Cosolvent_Head", r_n] = (
            cosolvent_data[f"{r_n}_Cosolvent_Head_integra_neg"][0]
            + cosolvent_data[f"{r_n}_Cosolvent_Head_integra_pos"][0]
        ) / 2
        output_df.loc["Cosolvent_Bulk", r_n] = (
            cosolvent_data[f"{r_n}_Cosolvent_Bulk_integra_neg"][0]
            + cosolvent_data[f"{r_n}_Cosolvent_Bulk_integra_pos"][0]
        ) / 2
        # output_df.loc["Cosolvent_Tail_len", r_n] = (
        #     cosolvent_data[f"{r_n}_Cosolvent_Tail_len_neg"][0]
        #     + cosolvent_data[f"{r_n}_Cosolvent_Tail_len_pos"][0]
        # ) / 2
        # output_df.loc["Cosolvent_Head_len", r_n] = (
        #     cosolvent_data[f"{r_n}_Cosolvent_Head_len_neg"][0]
        #     + cosolvent_data[f"{r_n}_Cosolvent_Head_len_pos"][0]
        # ) / 2
        output_df.loc["Bulk_Length", r_n] = (
            cosolvent_data[f"{r_n}_Cosolvent_Bulk_len_neg"][0]
            + cosolvent_data[f"{r_n}_Cosolvent_Bulk_len_pos"][0]
        ) / 2
        output_df.loc["Water in Lipid Head", r_n] = (
            cosolvent_data[f"{r_n}_water_head_integra_neg"][0]
            + cosolvent_data[f"{r_n}_water_head_integra_pos"][0]
        ) / 2
    output_df.loc["Dc"] = output_df.loc["C21C31"]
    output_df.loc["Dh"] = output_df.loc["90%Water"] - output_df.loc["C21C31"]
    output_df.loc["DB"] = 2 * output_df.loc["Dc"] + output_df.loc["Dh"]
    output_df.loc["Dpp"] = 2 * output_df.loc["P"]
    output_df.loc["Kp"] = (
        (output_df.loc["Cosolvent_Tail"] + output_df.loc["Cosolvent_Head"])
        / (output_df.loc["Dc"] + output_df.loc["Dh"])
    ) / (output_df.loc["Cosolvent_Bulk"] / output_df.loc["Bulk_Length"])
    output_df.loc["Ps"] = output_df.loc["Cosolvent_Tail"] / (
        output_df.loc["Cosolvent_Head"] + output_df.loc["Cosolvent_Tail"]
    )
    output_df.loc["Nw"] = (
        output_df.loc["Water in Lipid Head"]
        * output_df.loc["Average Box Area"]
        / lipid_perleaflet
    )
    return output_df
