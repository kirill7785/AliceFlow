object FormSolidLibMat: TFormSolidLibMat
  Left = 674
  Top = 261
  Caption = 'Solid Library Material'
  ClientHeight = 286
  ClientWidth = 498
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GBSelectMat: TGroupBox
    Left = 8
    Top = 16
    Width = 481
    Height = 265
    Caption = 'Please Select anyone SOLID Material'
    Color = clMoneyGreen
    ParentColor = False
    TabOrder = 0
    object LMessage1: TLabel
      Left = 16
      Top = 24
      Width = 316
      Height = 13
      Caption = 
        'Choose from realistic material properties depend on the temperat' +
        'ure'
    end
    object lblConduct: TLabel
      Left = 32
      Top = 88
      Width = 58
      Height = 13
      Caption = 'Conductivity'
    end
    object lblcapacity: TLabel
      Left = 32
      Top = 120
      Width = 41
      Height = 13
      Caption = 'Capacity'
    end
    object lblDensity: TLabel
      Left = 32
      Top = 160
      Width = 35
      Height = 13
      Caption = 'Density'
    end
    object lblvalDensity: TLabel
      Left = 112
      Top = 160
      Width = 3
      Height = 13
    end
    object lblDSI: TLabel
      Left = 168
      Top = 160
      Width = 34
      Height = 13
      Caption = 'kg/m!3'
    end
    object BApply: TButton
      Left = 368
      Top = 232
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 0
      OnClick = BApplyClick
    end
    object cbbSolidMatLib: TComboBox
      Left = 24
      Top = 48
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 1
      Text = 'Alumina Al2O3 polycrystalline aluminum oxide'
      OnChange = cbbSolidMatLibChange
      Items.Strings = (
        'Alumina Al2O3 polycrystalline aluminum oxide'
        'Si'
        'GaAs'
        'GaN'
        'SiC4H'
        'Sapphire'
        'Diamond'
        'MD40'
        'Au'
        'SiO2'
        'Cu'
        'Kovar'
        'BrassLS59_1_L'
        'Al_Duralumin'
        'AlN')
    end
    object btnviewL: TButton
      Left = 112
      Top = 88
      Width = 75
      Height = 25
      Caption = 'View'
      TabOrder = 2
      OnClick = btnviewLClick
    end
    object btncapacity: TButton
      Left = 112
      Top = 120
      Width = 75
      Height = 25
      Caption = 'View'
      TabOrder = 3
      OnClick = btncapacityClick
    end
  end
  object GroupBoxOrthotropy: TGroupBox
    Left = 232
    Top = 64
    Width = 185
    Height = 125
    Caption = 'Conductivity orthotropy multiplyer'
    TabOrder = 1
    object Labelx: TLabel
      Left = 16
      Top = 32
      Width = 7
      Height = 13
      Caption = 'X'
    end
    object Labely: TLabel
      Left = 16
      Top = 64
      Width = 7
      Height = 13
      Caption = 'Y'
    end
    object Labelz: TLabel
      Left = 16
      Top = 88
      Width = 7
      Height = 13
      Caption = 'Z'
    end
    object Editx: TEdit
      Left = 48
      Top = 24
      Width = 121
      Height = 21
      TabOrder = 0
    end
    object Edity: TEdit
      Left = 48
      Top = 51
      Width = 121
      Height = 21
      TabOrder = 1
    end
    object Editz: TEdit
      Left = 48
      Top = 78
      Width = 121
      Height = 21
      TabOrder = 2
    end
  end
end
