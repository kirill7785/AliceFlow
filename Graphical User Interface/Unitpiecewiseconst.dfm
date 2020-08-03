object Formpiecewiseconstant: TFormpiecewiseconstant
  Left = 0
  Top = 0
  Caption = 'Piecewise coonstant power law'
  ClientHeight = 281
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
  object Label2: TLabel
    Left = 119
    Top = 8
    Width = 64
    Height = 13
    Caption = 'time step (s);'
  end
  object Label3: TLabel
    Left = 216
    Top = 8
    Width = 83
    Height = 13
    Caption = 'power multiplyer.'
  end
  object Memopiecewiseconst: TMemo
    Left = 0
    Top = 32
    Width = 410
    Height = 241
    ScrollBars = ssVertical
    TabOrder = 0
  end
  object ComboBoxpiecewiseconst: TComboBox
    Left = 8
    Top = 8
    Width = 105
    Height = 21
    ItemIndex = 0
    TabOrder = 1
    Text = 'time (s);'
    Items.Strings = (
      'time (s);'
      'time duration (s);')
  end
end
