object FormFluidLibMat: TFormFluidLibMat
  Left = 370
  Top = 364
  AutoSize = True
  Caption = 'Fluid Library Material'
  ClientHeight = 265
  ClientWidth = 497
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnClose = FormClose
  PixelsPerInch = 96
  TextHeight = 13
  object GBFluidMatLib: TGroupBox
    Left = 0
    Top = 0
    Width = 497
    Height = 265
    Caption = 'Please select anyone FLUID material'
    Color = clMoneyGreen
    ParentColor = False
    TabOrder = 0
    object LFluidMatLib: TLabel
      Left = 16
      Top = 24
      Width = 316
      Height = 13
      Caption = 
        'Choose from realistic material properties depend on the temperat' +
        'ure'
    end
    object Label1: TLabel
      Left = 24
      Top = 80
      Width = 98
      Height = 13
      Caption = 'Thermal conductivity'
    end
    object Label2: TLabel
      Left = 209
      Top = 84
      Width = 42
      Height = 13
      Caption = 'W/(mxK)'
    end
    object Label3: TLabel
      Left = 24
      Top = 112
      Width = 67
      Height = 13
      Caption = 'Heat Capacity'
    end
    object Label4: TLabel
      Left = 209
      Top = 112
      Width = 40
      Height = 13
      Caption = 'J/(kgxK)'
    end
    object Label5: TLabel
      Left = 24
      Top = 144
      Width = 35
      Height = 13
      Caption = 'Density'
    end
    object LabelDensityValue: TLabel
      Left = 113
      Top = 144
      Width = 88
      Height = 13
      Caption = 'LabelDensityValue'
    end
    object Label6: TLabel
      Left = 219
      Top = 144
      Width = 34
      Height = 13
      Caption = 'kg/m!3'
    end
    object Label7: TLabel
      Left = 24
      Top = 176
      Width = 85
      Height = 13
      Caption = 'Dynamic Viscosity'
    end
    object Label8: TLabel
      Left = 209
      Top = 176
      Width = 23
      Height = 13
      Caption = 'Paxs'
    end
    object Label9: TLabel
      Left = 24
      Top = 216
      Width = 132
      Height = 13
      Caption = 'Linear expansion coefficient'
    end
    object Label10: TLabel
      Left = 255
      Top = 216
      Width = 18
      Height = 13
      Caption = '1/K'
    end
    object BApply: TButton
      Left = 392
      Top = 231
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 0
      OnClick = BApplyClick
    end
    object ComboBoxFluidLibMaterial: TComboBox
      Left = 16
      Top = 43
      Width = 185
      Height = 21
      ItemIndex = 0
      TabOrder = 1
      Text = 'Dry Air'
      OnChange = ComboBoxFluidLibMaterialChange
      Items.Strings = (
        'Dry Air'
        'Water Liquid')
    end
    object Buttonconductivity: TButton
      Left = 128
      Top = 72
      Width = 75
      Height = 25
      Caption = 'View'
      TabOrder = 2
      OnClick = ButtonconductivityClick
    end
    object ButtonHeatCapacity: TButton
      Left = 128
      Top = 103
      Width = 75
      Height = 25
      Caption = 'View'
      TabOrder = 3
      OnClick = ButtonHeatCapacityClick
    end
    object ButtonDynamicViscosity: TButton
      Left = 128
      Top = 168
      Width = 75
      Height = 25
      Caption = 'View'
      TabOrder = 4
      OnClick = ButtonDynamicViscosityClick
    end
    object Buttonbeta: TButton
      Left = 174
      Top = 211
      Width = 75
      Height = 25
      Caption = 'View'
      TabOrder = 5
      OnClick = ButtonbetaClick
    end
  end
end
