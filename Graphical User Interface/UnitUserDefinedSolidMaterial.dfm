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
      Height = 289
      Caption = 'Thermal-Stress'
      TabOrder = 4
      object GroupBoxPoissonRatio: TGroupBox
        Left = 16
        Top = 16
        Width = 289
        Height = 57
        Caption = 'Poisson Ratio'
        TabOrder = 0
        object EditPoissonRatio: TEdit
          Left = 16
          Top = 24
          Width = 121
          Height = 21
          TabOrder = 0
        end
      end
      object GroupBoxYoungModule: TGroupBox
        Left = 16
        Top = 79
        Width = 289
        Height = 66
        Caption = 'Young Module'
        TabOrder = 1
        object LabelYoungModule: TLabel
          Left = 143
          Top = 32
          Width = 21
          Height = 13
          Caption = 'GPa'
        end
        object EditYoungModule: TEdit
          Left = 16
          Top = 24
          Width = 121
          Height = 21
          TabOrder = 0
        end
      end
      object GroupBoxLinearExpansionCoefficient: TGroupBox
        Left = 16
        Top = 151
        Width = 289
        Height = 99
        Caption = 'Linear expansion coefficient'
        TabOrder = 2
        object LabelLinearExpansionCoefficient: TLabel
          Left = 143
          Top = 24
          Width = 47
          Height = 13
          Caption = '*1E-6 1/K'
        end
        object EditLinearExpansionKoefficient: TEdit
          Left = 16
          Top = 24
          Width = 121
          Height = 21
          TabOrder = 0
        end
      end
    end
  end
end
