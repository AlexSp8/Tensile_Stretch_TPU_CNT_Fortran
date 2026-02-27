import pandas as pd
from matplotlib import pyplot as plt
import os


working_dir = os.getcwd()


class Plot(object):

    def __init__(self, filenames: list, title: str, plot_directory_name="plots", plot_filename="plot_data"):

        self.filenames = filenames
        self.title = title
        self.dir_name = plot_directory_name
        self.plt_name = plot_filename
        self.plot_dir = os.path.join(working_dir, self.dir_name)

    ###############################################################
    def make_plot_directory(self):

        if not os.path.exists(self.plot_dir):
            os.mkdir(
                self.plot_dir
            )
            print("Folder generated")

        return None

    ###############################################################
    def plot(self):
        self.make_plot_directory()
        self.set_general_plot_confs()

        title = self.title
        plot_filename = self.plt_name

        for filename in self.filenames:

            if filename.__contains__("Stress"):
                self.plot_stress_strain(filename)
            elif filename.__contains__("size"):
                self.plot_specimen_size(filename)
            else:
                print(f"No plot function found for filename: {filename}")

        self.finalize()

    ###############################################################
    def finalize(self):
        # plt.show()
        plt.savefig(f"{self.plot_dir}/{self.plt_name}.png", dpi=300)
        # plt.savefig(f"{self.plot_dir}/{self.plt_name}.svg", dpi=1000, format="svg")
        plt.close()

    ###############################################################
    def set_general_plot_confs(self):
        plt.title(self.title)

    ###############################################################
    def plot_stress_strain(self, filename):

        if filename.__contains__("data/"):
            col_names = ['t', 'e_eng', 's_eng']
            df = pd.read_csv(filename, delim_whitespace=True, names=col_names, header=None)
            plt.scatter(df['e_eng'], df['s_eng'], label="Experiment", linewidth=1)
        elif filename.__contains__(".DAT"):
            col_names = ['t', 'e_eng', 'e', 's_eng', 's']
            df = pd.read_csv(filename, delim_whitespace=True, names=col_names, header=None)
            plt.plot(df['e_eng'], df['s_eng'], color='r', label="Numerical", linewidth=1)
            # plt.plot(df['e'], df['s'], color='r', label="Numerical", linewidth=1)
        else:
            print("Wrong data for Stress strain")

        # plt.yscale("log")
        # plt.xscale("log")
        plt.legend()

    ###############################################################
    def plot_specimen_size(self, filename):

        if filename.__contains__(".DAT"):
            col_names = ['t', 'xmin', 'z_xmin', 'ymin', 'z_ymin', 'zmax']
            df = pd.read_csv(filename, delim_whitespace=True, names=col_names, header=None)
            plt.plot(df['t'], df['xmin'], color='r', label="Xmin", linewidth=1)
            plt.plot(df['t'], df['ymin'], color='b', label="Ymin", linewidth=1)
        else:
            print("Wrong data for Stress strain")

        # plt.yscale("log")
        # plt.xscale("log")
        plt.legend()

data_folder = 'data/'

#Neat, CNT, CB, ABS
Specimen = 'CNT'

###############################################################
#Stress-strain
GD_i = []
GD_f = []

for i in range(0,1):

    files = []
    if Specimen == 'Neat':
        # file = data_folder+"Stress_neat_black.DTA"
        # files.append(file)
        # file = data_folder+"Stress_neat_blue.DTA"
        # files.append(file)
        file = data_folder+"Stress_neat.DTA"
        files.append(file)
    elif Specimen == 'CNT':
        # file = data_folder+"Stress_CNT_blue.DTA"
        # files.append(file)
        # file = data_folder+"Stress_CNT_red.DTA"
        # files.append(file)
        # file = data_folder+"Stress_CNT_black.DTA"
        # files.append(file)
        file = data_folder+"Stress_CNT2.DTA"
        files.append(file)
        # file = data_folder+"Stress_CNT3.DTA"
        # files.append(file)
        # file = data_folder+"Stress_CNT4.DTA"
        # files.append(file)
    elif Specimen == 'CB':
        file = data_folder+"Stress_CB_blue.DTA"
        files.append(file)
        file = data_folder+"Stress_CB_red.DTA"
        files.append(file)
        file = data_folder+"Stress_CB_black.DTA"
        files.append(file)
    elif Specimen == 'ABS':
        # file = data_folder+"Stress_ABS_hot.DTA"
        # files.append(file)
        # file = data_folder+"Stress_ABS_cold.DTA"
        # files.append(file)
        file = data_folder+"Stress_ABS_room.DTA"
        files.append(file)
    else:
        print('Wrong Specimen data selection!')
        
    file = "out/Stress.DAT"
    files.append(file)

    plotter = Plot(
        # filenames=[file_1.format(GD_i1 = GD_i[i], GD_f1 = GD_f[i]), \
        #             file_2.format(GD_i1 = GD_i[i], GD_f1 = GD_f[i])],
        filenames=files,
        plot_directory_name='Graphs',
        plot_filename="Stress-strain",
        title="Stress-strain"
    )
    plotter.plot()
    del plotter

###############################################################
#Size
for i in range(0,1):

    files = []
    file = "out/Specimen_size.DAT"
    files.append(file)

    plotter = Plot(
        filenames=files,
        plot_directory_name='Graphs',
        plot_filename="Specimen_size",
        title="Specimen size"
    )
    plotter.plot()
    del plotter
