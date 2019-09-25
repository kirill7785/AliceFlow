object FormTransientPowerSetting: TFormTransientPowerSetting
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Transient  Power  Setting in block'
  ClientHeight = 281
  ClientWidth = 171
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnClose = FormClose
  PixelsPerInch = 96
  TextHeight = 13
  object Panel1: TPanel
    Left = 0
    Top = 0
    Width = 171
    Height = 281
    TabOrder = 0
    object RadioGroupTimeDependPowerLow: TRadioGroup
      Left = 1
      Top = 1
      Width = 169
      Height = 137
      Caption = 'Time depent power law'
      ItemIndex = 0
      Items.Strings = (
        'it does not depend on time'
        'square wave'
        'square wave apparat'
        'hot cold double linear')
      TabOrder = 0
    end
    object Button1: TButton
      Left = 95
      Top = 248
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = Button1Click
    end
    object PanelTemperaturedependpower: TPanel
      Left = 10
      Top = 144
      Width = 161
      Height = 98
      TabOrder = 2
      object Label1: TLabel
        Left = 112
        Top = 40
        Width = 10
        Height = 13
        Caption = 'W'
      end
      object ComboBoxTemperaturedependpower: TComboBox
        Left = 8
        Top = 8
        Width = 145
        Height = 21
        ItemIndex = 0
        TabOrder = 0
        Text = 'Constant'
        OnChange = ComboBoxTemperaturedependpowerChange
        Items.Strings = (
          'Constant'
          'Piecewise linear')
      end
      object EditPower: TEdit
        Left = 8
        Top = 35
        Width = 89
        Height = 21
        TabOrder = 1
      end
      object Buttonpiecewisepower: TButton
        Left = 16
        Top = 64
        Width = 75
        Height = 25
        Caption = 'Edit'
        TabOrder = 2
        OnClick = ButtonpiecewisepowerClick
      end
    end
  end
end
