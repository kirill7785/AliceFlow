object FormOptimetric: TFormOptimetric
  Left = 0
  Top = 0
  Caption = 'Parameters and optimization'
  ClientHeight = 299
  ClientWidth = 635
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
    Left = 168
    Top = 8
    Width = 167
    Height = 13
    Caption = 'list of values '#8203#8203'separated by a space'
  end
  object EditListVariable: TEdit
    Left = 168
    Top = 32
    Width = 425
    Height = 21
    TabOrder = 0
  end
  object ComboBoxvar_id0: TComboBox
    Left = 17
    Top = 32
    Width = 145
    Height = 21
    TabOrder = 1
  end
  object ButtonRun: TButton
    Left = 518
    Top = 104
    Width = 75
    Height = 25
    Caption = 'Run'
    TabOrder = 2
    OnClick = ButtonRunClick
  end
  object ComboBoxvar_id1: TComboBox
    Left = 17
    Top = 59
    Width = 145
    Height = 21
    TabOrder = 3
  end
  object EditListVariable1: TEdit
    Left = 168
    Top = 59
    Width = 425
    Height = 21
    TabOrder = 4
  end
end
