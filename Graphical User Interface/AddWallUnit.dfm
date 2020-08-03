object AddWallForm: TAddWallForm
  Left = 208
  Top = 119
  AutoSize = True
  Caption = 'editing the properties of the wall'
  ClientHeight = 468
  ClientWidth = 425
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Panelglobalcontainer: TPanel
    Left = 0
    Top = 0
    Width = 425
    Height = 468
    Color = clMoneyGreen
    ParentBackground = False
    TabOrder = 0
    object RadioGroup1: TRadioGroup
      Left = 8
      Top = 11
      Width = 401
      Height = 62
      Caption = 'Size and properties wall'
      Columns = 3
      ItemIndex = 0
      Items.Strings = (
        'Info'
        'Geometry'
        'Properties')
      TabOrder = 0
      OnClick = RadioGroup1Click
    end
    object PanelInfo: TPanel
      Left = 8
      Top = 79
      Width = 401
      Height = 337
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 1
      object Lname: TLabel
        Left = 24
        Top = 24
        Width = 26
        Height = 13
        Caption = 'name'
      end
      object Ename: TEdit
        Left = 64
        Top = 21
        Width = 97
        Height = 21
        TabOrder = 0
      end
    end
    object Bapply: TButton
      Left = 309
      Top = 422
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 2
      OnClick = BapplyClick
    end
    object PanelGeometry: TPanel
      Left = 8
      Top = 79
      Width = 401
      Height = 337
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 3
      Visible = False
      object RadioGroupPlane: TRadioGroup
        Left = 16
        Top = 21
        Width = 209
        Height = 53
        Caption = 'Plane'
        Color = clMoneyGreen
        Columns = 3
        Items.Strings = (
          'XY'
          'XZ'
          'YZ')
        ParentBackground = False
        ParentColor = False
        TabOrder = 0
        OnClick = RadioGroupPlaneClick
      end
      object GroupBoxSize: TGroupBox
        Left = 8
        Top = 80
        Width = 305
        Height = 121
        Caption = 'Size'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 1
        object LxS: TLabel
          Left = 16
          Top = 24
          Width = 12
          Height = 13
          Caption = 'xS'
        end
        object LyS: TLabel
          Left = 16
          Top = 56
          Width = 12
          Height = 13
          Caption = 'yS'
        end
        object LzS: TLabel
          Left = 16
          Top = 88
          Width = 12
          Height = 13
          Caption = 'zS'
        end
        object LxE: TLabel
          Left = 158
          Top = 24
          Width = 12
          Height = 13
          Caption = 'xE'
        end
        object LyE: TLabel
          Left = 158
          Top = 56
          Width = 12
          Height = 13
          Caption = 'yE'
        end
        object LzE: TLabel
          Left = 158
          Top = 88
          Width = 12
          Height = 13
          Caption = 'zE'
        end
        object Label6: TLabel
          Left = 104
          Top = 24
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label7: TLabel
          Left = 104
          Top = 56
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label8: TLabel
          Left = 103
          Top = 88
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label9: TLabel
          Left = 255
          Top = 19
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label10: TLabel
          Left = 256
          Top = 56
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Label11: TLabel
          Left = 256
          Top = 88
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object ExS: TEdit
          Left = 40
          Top = 16
          Width = 57
          Height = 21
          TabOrder = 0
        end
        object EyS: TEdit
          Left = 40
          Top = 48
          Width = 57
          Height = 21
          TabOrder = 1
        end
        object EzS: TEdit
          Left = 40
          Top = 80
          Width = 57
          Height = 21
          TabOrder = 2
        end
        object ExE: TEdit
          Left = 176
          Top = 16
          Width = 73
          Height = 21
          TabOrder = 3
        end
        object EyE: TEdit
          Left = 176
          Top = 53
          Width = 73
          Height = 21
          TabOrder = 4
        end
        object EzE: TEdit
          Left = 176
          Top = 80
          Width = 73
          Height = 21
          TabOrder = 5
        end
      end
    end
    object PanelProperties: TPanel
      Left = 8
      Top = 80
      Width = 401
      Height = 336
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 4
      Visible = False
      object GroupBoxFLOW: TGroupBox
        Left = 7
        Top = 21
        Width = 153
        Height = 258
        Caption = 'FLOW Boundary Condition'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 0
        object RadioGroupflowtype: TRadioGroup
          Left = 8
          Top = 16
          Width = 121
          Height = 121
          Caption = 'Boundary Type'
          Color = clMoneyGreen
          ItemIndex = 0
          Items.Strings = (
            'specified velocity'
            'specified pressure'
            'symmetry'
            'opening')
          ParentBackground = False
          ParentColor = False
          TabOrder = 0
          OnClick = RadioGroupflowtypeClick
        end
        object GroupBoxvelcomp: TGroupBox
          Left = 3
          Top = 151
          Width = 137
          Height = 98
          Caption = 'Velocity component'
          Color = clMoneyGreen
          ParentBackground = False
          ParentColor = False
          TabOrder = 1
          object LabelVx: TLabel
            Left = 8
            Top = 20
            Width = 12
            Height = 13
            Caption = 'Vx'
          end
          object LabelVy: TLabel
            Left = 8
            Top = 48
            Width = 12
            Height = 13
            Caption = 'Vy'
          end
          object LabelVz: TLabel
            Left = 8
            Top = 72
            Width = 12
            Height = 13
            Caption = 'Vz'
          end
          object Label2: TLabel
            Left = 96
            Top = 64
            Width = 18
            Height = 13
            Caption = 'm/s'
          end
          object Label3: TLabel
            Left = 95
            Top = 45
            Width = 18
            Height = 13
            Caption = 'm/s'
          end
          object Label4: TLabel
            Left = 96
            Top = 24
            Width = 18
            Height = 13
            Caption = 'm/s'
          end
          object EditVx: TEdit
            Left = 32
            Top = 16
            Width = 57
            Height = 21
            TabOrder = 0
          end
          object EditVy: TEdit
            Left = 32
            Top = 40
            Width = 57
            Height = 21
            TabOrder = 1
          end
          object EditVz: TEdit
            Left = 32
            Top = 64
            Width = 57
            Height = 21
            TabOrder = 2
          end
        end
        object GroupBoxpressure: TGroupBox
          Left = 3
          Top = 163
          Width = 137
          Height = 81
          Caption = 'specified pressure'
          Color = clMoneyGreen
          ParentBackground = False
          ParentColor = False
          TabOrder = 2
          object Labelpressure: TLabel
            Left = 8
            Top = 24
            Width = 40
            Height = 13
            Caption = 'pressure'
          end
          object Label1: TLabel
            Left = 95
            Top = 48
            Width = 13
            Height = 13
            Caption = 'Pa'
          end
          object Editpress: TEdit
            Left = 8
            Top = 48
            Width = 81
            Height = 21
            TabOrder = 0
          end
        end
      end
      object GroupBoxtemper: TGroupBox
        Left = 167
        Top = 8
        Width = 217
        Height = 274
        Caption = 'Temperature Boundary Condition'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 1
        object RadioGroupBonConTemp: TRadioGroup
          Left = 8
          Top = 16
          Width = 201
          Height = 113
          Caption = 'Boundary Type'
          Color = clMoneyGreen
          ItemIndex = 0
          Items.Strings = (
            'Temperature'
            'homogeneous Neumann conditions'
            'Newton-Richman'
            'Stefan-Bolcman')
          ParentBackground = False
          ParentColor = False
          TabOrder = 0
          OnClick = RadioGroupBonConTempClick
        end
        object PaneltemperatureBC: TPanel
          Left = 12
          Top = 228
          Width = 192
          Height = 43
          Color = clMoneyGreen
          ParentBackground = False
          TabOrder = 1
          object LTemp: TLabel
            Left = 11
            Top = 17
            Width = 56
            Height = 13
            Caption = 'temperature'
          end
          object Label5: TLabel
            Left = 160
            Top = 16
            Width = 14
            Height = 13
            Caption = ' '#176#1057
          end
          object Etemp: TEdit
            Left = 73
            Top = 13
            Width = 81
            Height = 21
            TabOrder = 0
          end
        end
        object Panelemissivity: TPanel
          Left = 12
          Top = 135
          Width = 197
          Height = 82
          Color = clMoneyGreen
          ParentBackground = False
          TabOrder = 2
          Visible = False
          object Label12: TLabel
            Left = 5
            Top = 16
            Width = 44
            Height = 13
            Caption = 'emissivity'
          end
          object Label13: TLabel
            Left = 32
            Top = 43
            Width = 56
            Height = 13
            Caption = 'View Factor'
          end
          object Editemissivity: TEdit
            Left = 128
            Top = 16
            Width = 49
            Height = 21
            TabOrder = 0
            Text = '0.8'
          end
          object EditViewFactor: TEdit
            Left = 128
            Top = 43
            Width = 49
            Height = 21
            TabOrder = 1
            Text = '1.0'
          end
        end
      end
      object GroupBoxThermalStress: TGroupBox
        Left = 8
        Top = 285
        Width = 377
        Height = 45
        Caption = 'Thermal-Stress'
        TabOrder = 2
        object Labeldeformation: TLabel
          Left = 8
          Top = 24
          Width = 150
          Height = 13
          Caption = 'Deformation boundary condition'
        end
        object LabelForce: TLabel
          Left = 328
          Top = 24
          Width = 37
          Height = 13
          Caption = 'Newton'
        end
        object ComboBoxDeformationBoundaryConditon: TComboBox
          Left = 176
          Top = 16
          Width = 92
          Height = 21
          ItemIndex = 0
          TabOrder = 0
          Text = 'FREE'
          OnChange = ComboBoxDeformationBoundaryConditonChange
          Items.Strings = (
            'FREE'
            'X FIXIT'
            'Y FIXIT'
            'Z FIXIT'
            'XY FIXIT'
            'XZ FIXIT'
            'YZ FIXIT'
            'ALL FIXIT'
            'FORCE X'
            'FORCE Y'
            'FORCE Z')
        end
        object EditForce: TEdit
          Left = 274
          Top = 21
          Width = 47
          Height = 21
          TabOrder = 1
        end
      end
    end
  end
end
