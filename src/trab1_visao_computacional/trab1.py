import sys
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout, QLabel, QWidget, QLineEdit, QHBoxLayout, QVBoxLayout, QPushButton, QGroupBox
from PyQt5.QtGui import QDoubleValidator
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import pi, cos, sin


# Crie suas funções de translação, rotação, criação de referenciais, plotagem de setas e qualquer outra função que precisar

'''
COMECEI AQUI ######################################
'''


def move(dx, dy, dz):
    T = np.eye(4)
    T[0, -1] = dx
    T[1, -1] = dy
    T[2, -1] = dz
    return T


def z_rotation(graus):
    graus = float(graus)
    angle_rad = graus*pi/180
    rotation_matrix = np.array([[cos(angle_rad), -sin(angle_rad), 0, 0],
                               [sin(angle_rad), cos(angle_rad), 0, 0],
                               [0, 0, 1, 0],
                               [0, 0, 0, 1]])
    return rotation_matrix


def x_rotation(graus):
    graus = float(graus)
    angle_rad = graus*pi/180
    rotation_matrix = np.array([[1, 0, 0, 0],
                                [0, cos(angle_rad), -sin(angle_rad), 0],
                                [0, sin(angle_rad), cos(angle_rad), 0],
                                [0, 0, 0, 1]])
    return rotation_matrix


def y_rotation(graus):
    graus = float(graus)
    angle_rad = graus*pi/180
    rotation_matrix = np.array([[cos(angle_rad), 0, sin(angle_rad), 0],
                                [0, 1, 0, 0],
                                [-sin(angle_rad), 0, cos(angle_rad), 0],
                                [0, 0, 0, 1]])
    return rotation_matrix


def set_plot(ax=None, figure=None, lim=[-2, 2]):
    if figure == None:
        figure = plt.figure(figsize=(8, 8))
    if ax == None:
        ax = plt.axes(projection='3d')

    ax.set_title("camera reference")
    ax.set_xlim(lim)
    ax.set_xlabel("x axis")
    ax.set_ylim(lim)
    ax.set_ylabel("y axis")
    ax.set_zlim(lim)
    ax.set_zlabel("z axis")
    return ax

# adding quivers to the plot


def draw_arrows(point, base, axis, length=1.5):
    # The object base is a matrix, where each column represents the vector
    # of one of the axis, written in homogeneous coordinates (ax,ay,az,0)

    # Plot vector of x-axis
    axis.quiver(point[0], point[1], point[2], base[0, 0], base[1, 0],
                base[2, 0], color='red', pivot='tail',  length=length)
    # Plot vector of y-axis
    axis.quiver(point[0], point[1], point[2], base[0, 1], base[1, 1],
                base[2, 1], color='green', pivot='tail',  length=length)
    # Plot vector of z-axis
    axis.quiver(point[0], point[1], point[2], base[0, 2], base[1, 2],
                base[2, 2], color='blue', pivot='tail',  length=length)

    return axis


def cam_origem():
    # base vector values
    e1 = np.array([[1], [0], [0], [0]])  # X
    e2 = np.array([[0], [1], [0], [0]])  # Y
    e3 = np.array([[0], [0], [1], [0]])  # Z
    base = np.hstack((e1, e2, e3))

    # origin point
    point = np.array([[0], [0], [0], [1]])

    cam = np.hstack((base, point))

    return cam


def casa():
    house_coords = np.array([
        [0, 0, 0],
        [0, -10.0000, 0],
        [0, -10.0000, 12.0000],
        [0, -10.4000, 11.5000],
        [0, -5.0000, 16.0000],
        [0, 0, 12.0000],
        [0, 0.5000, 11.4000],
        [0, 0, 12.0000],
        [0, 0, 0],
        [-12.0000, 0, 0],
        [-12.0000, -5.0000, 0],
        [-12.0000, -10.0000, 0],
        [0, -10.0000, 0],
        [0, -10.0000, 12.0000],
        [-12.0000, -10.0000, 12.0000],
        [-12.0000, 0, 12.0000],
        [0, 0, 12.0000],
        [0, -10.0000, 12.0000],
        [0, -10.5000, 11.4000],
        [-12.0000, -10.5000, 11.4000],
        [-12.0000, -10.0000, 12.0000],
        [-12.0000, -5.0000, 16.0000],
        [0, -5.0000, 16.0000],
        [0, 0.5000, 11.4000],
        [-12.0000, 0.5000, 11.4000],
        [-12.0000, 0, 12.0000],
        [-12.0000, -5.0000, 16.0000],
        [-12.0000, -10.0000, 12.0000],
        [-12.0000, -10.0000, 0],
        [-12.0000, -5.0000, 0],
        [-12.0000, 0, 0],
        [-12.0000, 0, 12.0000],
        [-12.0000, 0, 0]
    ])

    window_coords_back = np.array([
        [-12.0000, -8.0, 4.0],
        [-12.0000, -2.0, 4.0],
        [-12.0000, -2.0, 8.0],
        [-12.0000, -8.0, 8.0],
        [-12.0000, -8.0, 4.0]
    ])

    window_coords_side = np.array([
        [-8.0, 0, 4.0],
        [-4.0, 0, 4.0],
        [-4.0, 0, 8.0],
        [-8.0, 0, 8.0],
        [-8.0, 0, 4.0]
    ])

    # Concatenar as coordenadas da casa e das duas janelas
    house = np.concatenate(
        (house_coords, window_coords_back, window_coords_side), axis=0)

    # Transpor a matriz para que cada coluna seja um ponto (X, Y, Z)
    house = np.transpose(house)

    # Adicionar um vetor de uns para representar as coordenadas homogêneas
    house = np.vstack([house, np.ones(np.size(house, 1))])

    return house


def mover_cam_pro_inicio(cam):
    Rx = x_rotation(-90)
    Rz = z_rotation(90)
    T = move(15, -5, 6)
    M = T@Rz@Rx
    cam = M@cam
    return cam


'''
Daqui pra baixo ja existia #######################
'''


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # definindo as variaveis
        self.set_variables()
        # Ajustando a tela
        self.setWindowTitle("Grid Layout")
        self.setGeometry(100, 100, 1280, 720)
        self.setup_ui()

    def set_variables(self):
        self.objeto_original = casa()  # modificar
        self.objeto = self.objeto_original
        self.cam_original = cam_origem()  # modificar
        self.cam = mover_cam_pro_inicio(self.cam_original)  # modificar
        self.px_base = 1280  # modificar
        self.px_altura = 720  # modificar
        self.dist_foc = 50  # modificar
        self.stheta = 0  # modificar
        self.ox = self.px_base/2  # modificar
        self.oy = self.px_altura/2  # modificar
        self.ccd = [36, 24]  # modificar
        self.projection_matrix = np.eye(3, 4)  # modificar

    def setup_ui(self):
        # Criar o layout de grade
        grid_layout = QGridLayout()

        # Criar os widgets
        line_edit_widget1 = self.create_world_widget("Ref mundo")
        line_edit_widget2 = self.create_cam_widget("Ref camera")
        line_edit_widget3 = self.create_intrinsic_widget("params instr")

        self.canvas = self.create_matplotlib_canvas()

        # Adicionar os widgets ao layout de grade
        grid_layout.addWidget(line_edit_widget1, 0, 0)
        grid_layout.addWidget(line_edit_widget2, 0, 1)
        grid_layout.addWidget(line_edit_widget3, 0, 2)
        grid_layout.addWidget(self.canvas, 1, 0, 1, 3)

        # Criar um widget para agrupar o botão de reset
        reset_widget = QWidget()
        reset_layout = QHBoxLayout()
        reset_widget.setLayout(reset_layout)

        # Criar o botão de reset vermelho
        reset_button = QPushButton("Reset")
        # Define um tamanho fixo para o botão (largura: 50 pixels, altura: 30 pixels)
        reset_button.setFixedSize(50, 30)
        style_sheet = """
            QPushButton {
                color : white ;
                background: rgba(255, 127, 130,128);
                font: inherit;
                border-radius: 5px;
                line-height: 1;
            }
        """
        reset_button.setStyleSheet(style_sheet)
        reset_button.clicked.connect(self.reset_canvas)

        # Adicionar o botão de reset ao layout
        reset_layout.addWidget(reset_button)

        # Adicionar o widget de reset ao layout de grade
        grid_layout.addWidget(reset_widget, 2, 0, 1, 3)

        # Criar um widget central e definir o layout de grade como seu layout
        central_widget = QWidget()
        central_widget.setLayout(grid_layout)

        # Definir o widget central na janela principal
        self.setCentralWidget(central_widget)

    def create_intrinsic_widget(self, title):
        # Criar um widget para agrupar os QLineEdit
        line_edit_widget = QGroupBox(title)
        line_edit_layout = QVBoxLayout()
        line_edit_widget.setLayout(line_edit_layout)

        # Criar um layout de grade para dividir os QLineEdit em 3 colunas
        grid_layout = QGridLayout()

        line_edits = []
        labels = ['n_pixels_base:', 'n_pixels_altura:', 'ccd_x:', 'ccd_y:',
                  'dist_focal:', 'sθ:']  # Texto a ser exibido antes de cada QLineEdit

        # Adicionar widgets QLineEdit com caixa de texto ao layout de grade
        for i in range(1, 7):
            line_edit = QLineEdit()
            label = QLabel(labels[i-1])
            validator = QDoubleValidator()  # Validador numérico
            # Aplicar o validador ao QLineEdit
            line_edit.setValidator(validator)
            grid_layout.addWidget(label, (i-1)//2, 2*((i-1) % 2))
            grid_layout.addWidget(line_edit, (i-1)//2, 2*((i-1) % 2) + 1)
            line_edits.append(line_edit)

        # Criar o botão de atualização
        update_button = QPushButton("Atualizar")

        # Você deverá criar, no espaço reservado ao final, a função self.update_params_intrinsc ou outra que você queira
        # Conectar a função de atualização aos sinais de clique do botão
        update_button.clicked.connect(
            lambda: self.update_params_intrinsc(line_edits))

        # Adicionar os widgets ao layout do widget line_edit_widget
        line_edit_layout.addLayout(grid_layout)
        line_edit_layout.addWidget(update_button)

        # Retornar o widget e a lista de caixas de texto
        return line_edit_widget

    def create_world_widget(self, title):
        # Criar um widget para agrupar os QLineEdit
        line_edit_widget = QGroupBox(title)
        line_edit_layout = QVBoxLayout()
        line_edit_widget.setLayout(line_edit_layout)

        # Criar um layout de grade para dividir os QLineEdit em 3 colunas
        grid_layout = QGridLayout()

        line_edits = []
        # Texto a ser exibido antes de cada QLineEdit
        labels = ['X(move):', 'X(angle):', 'Y(move):',
                  'Y(angle):', 'Z(move):', 'Z(angle):']

        # Adicionar widgets QLineEdit com caixa de texto ao layout de grade
        for i in range(1, 7):
            line_edit = QLineEdit()
            label = QLabel(labels[i-1])
            validator = QDoubleValidator()  # Validador numérico
            # Aplicar o validador ao QLineEdit
            line_edit.setValidator(validator)
            grid_layout.addWidget(label, (i-1)//2, 2*((i-1) % 2))
            grid_layout.addWidget(line_edit, (i-1)//2, 2*((i-1) % 2) + 1)
            line_edits.append(line_edit)

        # Criar o botão de atualização
        update_button = QPushButton("Atualizar")

        # Você deverá criar, no espaço reservado ao final, a função self.update_world ou outra que você queira
        # Conectar a função de atualização aos sinais de clique do botão
        update_button.clicked.connect(lambda: self.update_world(line_edits))

        # Adicionar os widgets ao layout do widget line_edit_widget
        line_edit_layout.addLayout(grid_layout)
        line_edit_layout.addWidget(update_button)

        # Retornar o widget e a lista de caixas de texto
        return line_edit_widget

    def create_cam_widget(self, title):
        # Criar um widget para agrupar os QLineEdit
        line_edit_widget = QGroupBox(title)
        line_edit_layout = QVBoxLayout()
        line_edit_widget.setLayout(line_edit_layout)

        # Criar um layout de grade para dividir os QLineEdit em 3 colunas
        grid_layout = QGridLayout()

        line_edits = []
        # Texto a ser exibido antes de cada QLineEdit
        labels = ['X(move):', 'X(angle):', 'Y(move):',
                  'Y(angle):', 'Z(move):', 'Z(angle):']

        # Adicionar widgets QLineEdit com caixa de texto ao layout de grade
        for i in range(1, 7):
            line_edit = QLineEdit()
            label = QLabel(labels[i-1])
            validator = QDoubleValidator()  # Validador numérico
            # Aplicar o validador ao QLineEdit
            line_edit.setValidator(validator)
            grid_layout.addWidget(label, (i-1)//2, 2*((i-1) % 2))
            grid_layout.addWidget(line_edit, (i-1)//2, 2*((i-1) % 2) + 1)
            line_edits.append(line_edit)

        # Criar o botão de atualização
        update_button = QPushButton("Atualizar")

        # Você deverá criar, no espaço reservado ao final, a função self.update_cam ou outra que você queira
        # Conectar a função de atualização aos sinais de clique do botão
        update_button.clicked.connect(lambda: self.update_cam(line_edits))

        # Adicionar os widgets ao layout do widget line_edit_widget
        line_edit_layout.addLayout(grid_layout)
        line_edit_layout.addWidget(update_button)

        # Retornar o widget e a lista de caixas de texto
        return line_edit_widget

    def create_matplotlib_canvas(self):
        # Criar um widget para exibir os gráficos do Matplotlib
        canvas_widget = QWidget()
        canvas_layout = QHBoxLayout()
        canvas_widget.setLayout(canvas_layout)

        # Criar um objeto FigureCanvas para exibir o gráfico 2D
        self.fig1, self.ax1 = plt.subplots()
        self.ax1.set_title("Imagem")
        self.canvas1 = FigureCanvas(self.fig1)

        # Falta acertar os limites do eixo X
        self.ax1.set_xlim([0, self.px_base])

        # Falta acertar os limites do eixo Y
        self.ax1.set_ylim([self.px_altura, 0])

        # Você deverá criar a função de projeção
        object_2d = self.projection_2d()

        # Falta plotar o object_2d que retornou da projeção
        self.ax1.plot(object_2d[0, :], object_2d[1, :])

        self.ax1.grid('True')
        self.ax1.set_aspect('equal')
        canvas_layout.addWidget(self.canvas1)

        # Criar um objeto FigureCanvas para exibir o gráfico 3D
        self.fig2 = plt.figure()
        self.ax2 = self.fig2.add_subplot(111, projection='3d')

        # Falta plotar o seu objeto 3D e os referenciais da câmera e do mundo
        self.ax2 = set_plot(ax=self.ax2, lim=[-15, 30])
        self.ax2.plot3D(self.objeto[0, :],
                        self.objeto[1, :], self.objeto[2, :], 'red')
        draw_arrows(self.cam[:, -1], self.cam[:, 0:3], self.ax2)

        self.canvas2 = FigureCanvas(self.fig2)
        canvas_layout.addWidget(self.canvas2)

        # Retornar o widget de canvas
        return canvas_widget

    # Você deverá criar as suas funções aqui

    def update_params_intrinsc(self, line_edits):
        px_b = line_edits[0].text()
        px_alt = line_edits[1].text()
        ccd_x = line_edits[2].text()
        ccd_y = line_edits[3].text()
        dist_f = line_edits[4].text()
        s_theta = line_edits[5].text()

        if px_b:
            self.px_base = float(px_b)
        if px_alt:
            self.px_altura = float(px_alt)
        if ccd_x:
            self.ccd[0] = float(ccd_x)
        if ccd_y:
            self.ccd[1] = float(ccd_y)
        if dist_f:
            self.dist_foc = float(dist_f)
        if s_theta:
            self.stheta = float(s_theta)

        self.update_canvas()
        for line_edit in line_edits:
            line_edit.clear()

    def update_world(self, line_edits):
        x_deslocamento = line_edits[0].text()
        x_angulo = line_edits[1].text()
        y_deslocamento = line_edits[2].text()
        y_angulo = line_edits[3].text()
        z_deslocamento = line_edits[4].text()
        z_angulo = line_edits[5].text()

        if x_deslocamento:
            dx = x_deslocamento
        else:
            dx = 0

        if y_deslocamento:
            dy = y_deslocamento
        else:
            dy = 0

        if z_deslocamento:
            dz = z_deslocamento
        else:
            dz = 0

        translacao = move(dx, dy, dz)
        self.cam = translacao@self.cam

        if x_angulo:
            self.cam = x_rotation(x_angulo) @ self.cam

        if y_angulo:
            self.cam = y_rotation(y_angulo)@self.cam

        if z_angulo:
            self.cam = z_rotation(z_angulo) @ self.cam

        self.update_canvas()
        for line_edit in line_edits:
            line_edit.clear()

    def update_cam(self, line_edits):
        x_deslocamento = line_edits[0].text()
        x_angulo = line_edits[1].text()
        y_deslocamento = line_edits[2].text()
        y_angulo = line_edits[3].text()
        z_deslocamento = line_edits[4].text()
        z_angulo = line_edits[5].text()

        if x_deslocamento:
            dx = x_deslocamento
        else:
            dx = 0

        if y_deslocamento:
            dy = y_deslocamento
        else:
            dy = 0

        if z_deslocamento:
            dz = z_deslocamento
        else:
            dz = 0

        trans = move(dx, dy, dz)
        self.cam = self.cam@trans@self.cam_original

        if x_angulo:
            self.cam = self.cam@x_rotation(x_angulo)@self.cam_original

        if y_angulo:
            self.cam = self.cam@y_rotation(y_angulo)@self.cam_original

        if z_angulo:
            self.cam = self.cam@z_rotation(z_angulo)@self.cam_original

        self.update_canvas()
        for line_edit in line_edits:
            line_edit.clear()

    def projection_2d(self):
        G = np.linalg.inv(self.cam)
        K = self.generate_intrinsic_params_matrix()

        '''
              [x]   [F*sx F*s_theta Ox] [1 0 0 0] [R T]^-1
        lambda[y] = [0     F*sy     Oy] [0 1 0 0] [0 1]
              [1]   [0      0        1] [0 0 1 0]
        '''
        object_2d = K@self.projection_matrix@G@self.objeto
        object_2d[0, :] = object_2d[0, :]/object_2d[2, :]
        object_2d[1, :] = object_2d[1, :]/object_2d[2, :]
        object_2d[2, :] = object_2d[2, :]/object_2d[2, :]
        return object_2d

    def generate_intrinsic_params_matrix(self):
        f = self.dist_foc
        w_p = self.px_base
        w_mm = self.ccd[0]

        h_p = self.px_altura
        h_mm = self.ccd[1]

        sx = w_p/w_mm
        sy = h_p/h_mm
        s_theta = self.stheta
        Ox = w_p/2
        Oy = h_p/2

        K = np.array([[f*sx, f*s_theta, Ox],
                      [0, f*sy, Oy],
                      [0, 0, 1]])

        return K

    def update_canvas(self):
        plt.close('all')
        object_2d = self.projection_2d()
        self.ax1.clear()
        self.ax1.set_xlim([0, self.px_base])
        self.ax1.set_ylim([self.px_altura, 0])
        self.ax1.plot(object_2d[0, :], object_2d[1, :])
        self.ax1.grid(True)
        self.ax1.set_aspect('equal')

        self.ax2.clear()
        self.ax2 = set_plot(ax=self.ax2, lim=[-15, 20])
        # draw_arrows(point,base,ax0)
        draw_arrows(self.cam[:, -1], self.cam[:, 0:3], self.ax2)
        self.ax2.plot3D(self.objeto[0, :],
                        self.objeto[1, :],
                        self.objeto[2, :],
                        'red')

        self.canvas1.draw()
        self.canvas2.draw()
        self.canvas.layout().itemAt(1).widget().draw()

    def reset_canvas(self):
        self.set_variables()
        self.update_canvas()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
