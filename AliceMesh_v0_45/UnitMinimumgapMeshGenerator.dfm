object FormMinimum_gap: TFormMinimum_gap
  Left = 0
  Top = 0
  Caption = 'Minimum gap'
  ClientHeight = 171
  ClientWidth = 232
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
    Left = 16
    Top = 16
    Width = 157
    Height = 13
    Caption = 'Minimum gap for mesh generator'
  end
  object Label2: TLabel
    Left = 24
    Top = 48
    Width = 6
    Height = 13
    Caption = 'X'
  end
  object Label3: TLabel
    Left = 168
    Top = 48
    Width = 8
    Height = 13
    Caption = 'm'
  end
  object Label4: TLabel
    Left = 24
    Top = 80
    Width = 6
    Height = 13
    Caption = 'Y'
  end
  object Label5: TLabel
    Left = 24
    Top = 112
    Width = 6
    Height = 13
    Caption = 'Z'
  end
  object Label6: TLabel
    Left = 168
    Top = 75
    Width = 8
    Height = 13
    Caption = 'm'
  end
  object Label7: TLabel
    Left = 168
    Top = 102
    Width = 8
    Height = 13
    Caption = 'm'
  end
  object EditX: TEdit
    Left = 36
    Top = 45
    Width = 121
    Height = 21
    TabOrder = 0
  end
  object Button1: TButton
    Left = 136
    Top = 138
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 1
    OnClick = Button1Click
  end
  object EditY: TEdit
    Left = 36
    Top = 72
    Width = 121
    Height = 21
    TabOrder = 2
  end
  object EditZ: TEdit
    Left = 36
    Top = 99
    Width = 121
    Height = 21
    TabOrder = 3
  end
end
