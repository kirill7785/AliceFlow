object FormUserDefinedSolidMat: TFormUserDefinedSolidMat
  Left = 306
  Top = 152
  Caption = 'User-Defined Solid Material '
  ClientHeight = 567
  ClientWidth = 614
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object PanelSOLID: TPanel
    Left = 0
    Top = 0
    Width = 617
    Height = 569
    Color = clMoneyGreen
    TabOrder = 0
    object GBuserproperties: TGroupBox
      Left = 8
      Top = 136
      Width = 273
      Height = 385
      Caption = 'user solid properties'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 0
      object GroupBoxRho: TGroupBox
        Left = 8
        Top = 16
        Width = 257
        Height = 49
        Caption = 'Density (Rho)'
        TabOrder = 0
        object LSIrho: TLabel
          Left = 208
          Top = 16
          Width = 37
          Height = 13
          Caption = 'kg/m^3'
        end
        object Erho: TEdit
          Left = 8
          Top = 16
          Width = 193
          Height = 21
          TabOrder = 0
        end
      end
      object GroupBoxCp: TGroupBox
        Left = 8
        Top = 72
        Width = 257
        Height = 89
        Caption = 'Heat Capacity (Cp)'
        TabOrder = 1
        object LSICp: TLabel
          Left = 159
          Top = 56
          Width = 39
          Height = 13
          Caption = 'J/(kg*K)'
        end
        object ECp: TEdit
          Left = 8
          Top = 56
          Width = 145
          Height = 21
          TabOrder = 0
        end
        object ComboBoxheatcapacitytype: TComboBox
          Left = 8
          Top = 24
          Width = 102
          Height = 21
          ItemIndex = 0
          TabOrder = 1
          Text = 'Constant'
          OnChange = ComboBoxheatcapacitytypeChange
          Items.Strings = (
            'Constant'
            'Piecewise')
        end
        object Buttonheatcapacitypiecewise: TButton
          Left = 126
          Top = 25
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 2
          OnClick = ButtonheatcapacitypiecewiseClick
        end
      end
      object GroupBoxLam: TGroupBox
        Left = 13
        Top = 167
        Width = 252
        Height = 89
        Caption = 'Conductivity (Lam)'
        TabOrder = 2
        object LSILam: TLabel
          Left = 166
          Top = 56
          Width = 41
          Height = 13
          Caption = 'W/(m*K)'
        end
        object ELam: TEdit
          Left = 15
          Top = 51
          Width = 145
          Height = 21
          TabOrder = 0
        end
        object ComboBoxconductivitytype: TComboBox
          Left = 15
          Top = 24
          Width = 90
          Height = 21
          ItemIndex = 0
          TabOrder = 1
          Text = 'Constant'
          OnChange = ComboBoxconductivitytypeChange
          Items.Strings = (
            'Constant'
            'Piecewise')
        end
        object Buttonconductiviypiecewise: TButton
          Left = 126
          Top = 20
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 2
          OnClick = ButtonconductiviypiecewiseClick
        end
      end
      object GroupBoxOrthotropy: TGroupBox
        Left = 16
        Top = 262
        Width = 241
        Height = 114
        Caption = 'Orthotropy conductivity multiplyer'
        TabOrder = 3
        object Labelx: TLabel
          Left = 40
          Top = 32
          Width = 7
          Height = 13
          Caption = 'X'
        end
        object Labely: TLabel
          Left = 40
          Top = 59
          Width = 7
          Height = 13
          Caption = 'Y'
        end
        object Labelz: TLabel
          Left = 40
          Top = 88
          Width = 7
          Height = 13
          Caption = 'Z'
        end
        object Editmultx: TEdit
          Left = 72
          Top = 24
          Width = 121
          Height = 21
          Color = clWhite
          TabOrder = 0
        end
        object Editmulty: TEdit
          Left = 72
          Top = 51
          Width = 121
          Height = 21
          TabOrder = 1
        end
        object Editmultz: TEdit
          Left = 72
          Top = 78
          Width = 121
          Height = 21
          TabOrder = 2
        end
      end
    end
    object BApply: TButton
      Left = 206
      Top = 527
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = BApplyClick
    end
    object GroupBoxsmn: TGroupBox
      Left = 8
      Top = 72
      Width = 273
      Height = 57
      Caption = 'solid material name'
      TabOrder = 2
      object EMatName: TEdit
        Left = 8
        Top = 24
        Width = 257
        Height = 21
        TabOrder = 0
      end
    end
    object GroupBoxtippattern: TGroupBox
      Left = 8
      Top = 8
      Width = 273
      Height = 57
      Caption = 'tip pattern'
      TabOrder = 3
      object CBsolidmat: TComboBox
        Left = 8
        Top = 24
        Width = 257
        Height = 21
        TabOrder = 0
        OnChange = CBsolidmatChange
        Items.Strings = (
          'no pattern'
          'Al-Duralumin'
          'GaN'
          'SiC-4H'
          'SOLDER_Au80Sn20'
          'Cu'
          'MD-40'
          'GaAs'
          'Au'
          'SiO2'
          'Si'
          'kovar'
          'Alumina'
          'Brass'
          'Polyimide'
          'Sapphire'
          'Glue_ECHES'
          'BeO'
          'Ag'
          'Diamond'
          'Si3N4'
          'KPT8'
          'Polyurethane_foam'
          'SOLDER_PbSn2Ag2.5'
          'SOLDER_SnAg25Sb10'
          'SOLDER_SnPb36Ag2'
          'GLUE_Ablebond_3230'
          'GLUE_Ablebond_8290'
          'GLUE_CRM_1033B'
          'SOLDER_Au88Ge12'
          'kanifoul'
          'Polystyrene_rigid_R12'
          'FR4'
          'Polystyrene_Typical'
          'air_solid'
          'indium'
          'VK_87_Ceramics'
          'AlN'
          'vacuum'
          'rogers5870'
          'rogers6002'
          'rogers5880'
          'rohacellhf51D')
      end
    end
    object GroupBoxThermalStress: TGroupBox
      Left = 287
      Top = 8
      Width = 322
      Height = 544
      Caption = 'Thermal-Stress'
      TabOrder = 4
      object GroupBoxPoissonRatio: TGroupBox
        Left = 16
        Top = 16
        Width = 289
        Height = 137
        Caption = 'Poisson Ratio'
        TabOrder = 0
        object Label5: TLabel
          Left = 16
          Top = 83
          Width = 158
          Height = 13
          Caption = 'Orthotropy Poisson ratio multiplyer'
        end
        object Label6: TLabel
          Left = 16
          Top = 102
          Width = 10
          Height = 13
          Caption = 'xy'
        end
        object Label7: TLabel
          Left = 95
          Top = 105
          Width = 10
          Height = 13
          Caption = 'xz'
        end
        object Label8: TLabel
          Left = 174
          Top = 105
          Width = 10
          Height = 13
          Caption = 'yz'
        end
        object EditPoissonRatio: TEdit
          Left = 16
          Top = 56
          Width = 121
          Height = 21
          TabOrder = 0
        end
        object ComboBoxPoissonratio: TComboBox
          Left = 16
          Top = 20
          Width = 129
          Height = 21
          ItemIndex = 0
          TabOrder = 1
          Text = 'Constant'
          OnClick = ComboBoxPoissonratioClick
          Items.Strings = (
            'Constant'
            'Piecewise')
        end
        object ButtonPoissonratio: TButton
          Left = 151
          Top = 18
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 2
          OnClick = ButtonPoissonratioClick
        end
        object Editnuxy: TEdit
          Left = 32
          Top = 102
          Width = 57
          Height = 21
          TabOrder = 3
          Text = '1.0'
        end
        object Editnuxz: TEdit
          Left = 111
          Top = 102
          Width = 57
          Height = 21
          TabOrder = 4
          Text = '1.0'
        end
        object Editnuyz: TEdit
          Left = 190
          Top = 102
          Width = 60
          Height = 21
          TabOrder = 5
          Text = '1.0'
        end
      end
      object GroupBoxYoungModule: TGroupBox
        Left = 16
        Top = 187
        Width = 289
        Height = 131
        Caption = 'Young Module'
        TabOrder = 1
        object LabelYoungModule: TLabel
          Left = 143
          Top = 55
          Width = 21
          Height = 13
          Caption = 'GPa'
        end
        object Label12: TLabel
          Left = 16
          Top = 78
          Width = 167
          Height = 13
          Caption = 'Orthotropy Young Module multiplyer'
        end
        object Label13: TLabel
          Left = 16
          Top = 104
          Width = 7
          Height = 13
          Caption = 'X'
        end
        object Label14: TLabel
          Left = 96
          Top = 104
          Width = 7
          Height = 13
          Caption = 'Y'
        end
        object Label15: TLabel
          Left = 176
          Top = 104
          Width = 7
          Height = 13
          Caption = 'Z'
        end
        object EditYoungModule: TEdit
          Left = 16
          Top = 51
          Width = 121
          Height = 21
          TabOrder = 0
        end
        object ButtonYoungModule: TButton
          Left = 152
          Top = 24
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 1
          OnClick = ButtonYoungModuleClick
        end
        object ComboBoxYoungModule: TComboBox
          Left = 16
          Top = 24
          Width = 130
          Height = 21
          ItemIndex = 0
          TabOrder = 2
          Text = 'Constant'
          OnClick = ComboBoxYoungModuleClick
          Items.Strings = (
            'Constant'
            'Piecewise')
        end
        object EditEx: TEdit
          Left = 29
          Top = 97
          Width = 52
          Height = 21
          TabOrder = 3
          Text = '1.0'
        end
        object EditEy: TEdit
          Left = 109
          Top = 97
          Width = 52
          Height = 21
          TabOrder = 4
          Text = '1.0'
        end
        object EditEz: TEdit
          Left = 189
          Top = 96
          Width = 60
          Height = 21
          TabOrder = 5
          Text = '1.0'
        end
      end
      object GroupBoxLinearExpansionCoefficient: TGroupBox
        Left = 16
        Top = 411
        Width = 289
        Height = 130
        Caption = 'Linear expansion coefficient'
        TabOrder = 2
        object LabelLinearExpansionCoefficient: TLabel
          Left = 175
          Top = 56
          Width = 47
          Height = 13
          Caption = '*1E-6 1/K'
        end
        object Label1: TLabel
          Left = 24
          Top = 100
          Width = 7
          Height = 13
          Caption = 'X'
        end
        object Label2: TLabel
          Left = 100
          Top = 100
          Width = 7
          Height = 13
          Caption = 'Y'
        end
        object Label3: TLabel
          Left = 179
          Top = 100
          Width = 7
          Height = 13
          Caption = 'Z'
        end
        object Label4: TLabel
          Left = 16
          Top = 78
          Width = 226
          Height = 13
          Caption = 'Orthotropy linear expansion coefficient multiplyer'
        end
        object EditLinearExpansionKoefficient: TEdit
          Left = 16
          Top = 51
          Width = 153
          Height = 21
          TabOrder = 0
        end
        object ComboBoxlinearExpansion: TComboBox
          Left = 16
          Top = 24
          Width = 121
          Height = 21
          ItemIndex = 0
          TabOrder = 1
          Text = 'Constant'
          OnChange = ComboBoxlinearExpansionChange
          Items.Strings = (
            'Constant'
            'Piecewise')
        end
        object ButtonEditlinearexpansioncoefficient: TButton
          Left = 152
          Top = 16
          Width = 75
          Height = 25
          Caption = 'Edit'
          TabOrder = 2
          OnClick = ButtonEditlinearexpansioncoefficientClick
        end
        object EditbetaX: TEdit
          Left = 37
          Top = 97
          Width = 57
          Height = 21
          TabOrder = 3
          Text = '1.0'
        end
        object EditbetaY: TEdit
          Left = 113
          Top = 98
          Width = 60
          Height = 21
          TabOrder = 4
          Text = '1.0'
        end
        object EditbetaZ: TEdit
          Left = 203
          Top = 97
          Width = 60
          Height = 21
          TabOrder = 5
          Text = '1.0'
        end
      end
      object GroupBoxShearModulus: TGroupBox
        Left = 16
        Top = 324
        Width = 289
        Height = 81
        Caption = 'Shear Modulus GPa'
        TabOrder = 3
        object CheckBoxShearModulus: TCheckBox
          Left = 8
          Top = 16
          Width = 97
          Height = 17
          Caption = 'Active'
          TabOrder = 0
          OnClick = CheckBoxShearModulusClick
        end
        object PanelShearModulus: TPanel
          Left = 16
          Top = 39
          Width = 270
          Height = 36
          Color = clMoneyGreen
          ParentBackground = False
          TabOrder = 1
          object Label9: TLabel
            Left = 8
            Top = 16
            Width = 18
            Height = 13
            Caption = 'Gxy'
          end
          object Label10: TLabel
            Left = 96
            Top = 16
            Width = 18
            Height = 13
            Caption = 'Gyz'
          end
          object Label11: TLabel
            Left = 183
            Top = 16
            Width = 18
            Height = 13
            Caption = 'Gxz'
          end
          object EditGxy: TEdit
            Left = 32
            Top = 10
            Width = 57
            Height = 21
            TabOrder = 0
            Text = '1.0'
          end
          object EditGyz: TEdit
            Left = 120
            Top = 8
            Width = 57
            Height = 21
            TabOrder = 1
            Text = '1.0'
          end
          object EditGxz: TEdit
            Left = 207
            Top = 8
            Width = 54
            Height = 21
            TabOrder = 2
            Text = '1.0'
          end
        end
      end
    end
    object ButtonCancel: TButton
      Left = 106
      Top = 527
      Width = 75
      Height = 25
      Caption = 'Cancel'
      TabOrder = 5
      OnClick = ButtonCancelClick
    end
  end
end
