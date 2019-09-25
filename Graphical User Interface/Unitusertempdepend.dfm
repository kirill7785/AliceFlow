object Formusertempdepend: TFormusertempdepend
  Left = 0
  Top = 0
  Caption = 'Formusertempdepend'
  ClientHeight = 291
  ClientWidth = 418
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Label1: TLabel
    Left = 8
    Top = 8
    Width = 31
    Height = 13
    Caption = 'Label1'
  end
  object Label2: TLabel
    Left = 8
    Top = 32
    Width = 31
    Height = 13
    Caption = 'Label2'
  end
  object Label3: TLabel
    Left = 8
    Top = 56
    Width = 31
    Height = 13
    Caption = 'Label3'
  end
  object Label4: TLabel
    Left = 8
    Top = 80
    Width = 31
    Height = 13
    Caption = 'Label4'
  end
  object Label5: TLabel
    Left = 216
    Top = 215
    Width = 31
    Height = 13
    Caption = 'Label5'
  end
  object Label6: TLabel
    Left = 16
    Top = 224
    Width = 91
    Height = 13
    Caption = 'Temperature units '
  end
  object Memopiecewiseproperties: TMemo
    Left = 8
    Top = 104
    Width = 402
    Height = 105
    ScrollBars = ssVertical
    TabOrder = 0
  end
  object ButtonApply: TButton
    Left = 304
    Top = 232
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 1
    OnClick = ButtonApplyClick
  end
  object ButtonView: TButton
    Left = 200
    Top = 234
    Width = 75
    Height = 25
    Caption = 'View'
    TabOrder = 2
    OnClick = ButtonViewClick
  end
  object ComboBoxtemperatureUnit: TComboBox
    Left = 112
    Top = 224
    Width = 41
    Height = 21
    ItemIndex = 0
    TabOrder = 3
    Text = ' '#176#1057
    Items.Strings = (
      ' '#176#1057
      'K')
  end
end
