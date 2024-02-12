#Visualize


import pandas as pd
import numpy as np
import os
import PySimpleGUI as sg
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.use('TkAgg')
brend_df = pd.read_csv("CSV/brend.csv",encoding ="shift-jis")

properties = brend_df.columns.tolist()[1:] #官能基の情報のみ

# GUIのレイアウト
sg.theme("DefaultNoMOreNagging")

layout = [[sg.Text('処方管理システム/Visualize', font=('Arial',20,"bold"))],
          [sg.Text("")],

          [sg.Text("X",size = (15,1)),
          sg.Combo(properties,
                   size=(30, 1),
                   key="-x-",
                   text_color='black',
                   background_color='honeydew'
                  )],

          [sg.Text("Y",size=(15, 1)),
          sg.Combo(properties,
                   size=(30, 1),
                   key="-y-",
                   text_color='black',
                   background_color='honeydew'
                  )],

         [sg.Text("Z",size=(15, 1)),
          sg.Combo(properties,
                   size=(30, 1),
                   key="-z-",
                   text_color='black',
                   background_color='honeydew'
                  )],

         [sg.Text("Marker",size=(15, 1)),
          sg.Combo(properties,
                   size=(30, 1),
                   key="marker",
                   text_color='black',
                   background_color='honeydew'
                  )],

          [sg.Push(),sg.Button("2次元散布図描画",size = (30,1))],
          [sg.Push(),sg.Button("3次元散布図描画",size = (30,1))],
          [sg.Push(),sg.Button("マルチプロット描画",size = (30,1))],
          


         ]



# ウインドウの出現位置を指定
win_location = (0, 0)

# モニターの解像度を取得
screen_width, screen_height = sg.Window.get_screen_size()

# ウィンドウのサイズをモニターの半分に設定
window_size = (int(screen_width*0.5),int(screen_height*0.5))

window = sg.Window("3dscatter",
                   layout,
                   size=window_size,
                   resizable=True,
                   location=win_location)


while True:

    event, values = window.read()

    x_name = values["-x-"]
    y_name = values["-y-"]
    z_name = values["-z-"]
    m_name = values["marker"]



    if event == sg.WINDOW_CLOSED:
        break
    if event == "2次元散布図描画":
        from modules import function as fc

        tag =brend_df["Name"].tolist()
        x = brend_df[x_name].tolist()
        y = brend_df[y_name].tolist()
        m = brend_df[m_name].tolist()

        labels = [x_name,y_name]
        fc.scatter2d(x,y,m,tag,labels)
        
    if event == "3次元散布図描画":
        from modules import function as fc

        tag =brend_df["Name"].tolist()
        x = brend_df[x_name].tolist()
        y = brend_df[y_name].tolist()
        z = brend_df[z_name].tolist()
        m = brend_df[m_name].tolist()

        labels = [x_name,y_name,z_name]
        fc.scatter3d(x,y,z,m,tag,labels)

    if event == "マルチプロット描画":
            
        num_rows = 5
        num_cols = 5
        
        # ループで10個のサブプロットを生成して表示
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 8))
        
        for ax, elem in zip(axes.flat,brend_df.columns[1:26]):
            x = brend_df[elem].tolist()
            m = brend_df[m_name].tolist()
            ax.scatter(x, m, s=3, marker='o', color='blue')  # 仮のデータをプロット
            ax.set_title(elem,fontsize=10)
            ax.set_ylabel('Target',fontsize=10)
        # サブプロット間のスペースを調整
        plt.tight_layout()
        plt.get_current_fig_manager().window.wm_geometry("+0+0")
        plt.show(block=False)


window.close()