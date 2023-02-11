object FormMoveSTLposition: TFormMoveSTLposition
  Left = 0
  Top = 0
  Caption = 'Move STL-position'
  ClientHeight = 152
  ClientWidth = 228
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
    Top = 8
    Width = 100
    Height = 13
    Caption = 'Move all STL-position'
  end
  object Label2: TLabel
    Left = 16
    Top = 40
    Width = 6
    Height = 13
    Caption = 'X'
  end
  object Label3: TLabel
    Left = 160
    Top = 40
    Width = 8
    Height = 13
    Caption = 'm'
  end
  object Label4: TLabel
    Left = 16
    Top = 72
    Width = 6
    Height = 13
    Caption = 'Y'
  end
  object Label5: TLabel
    Left = 160
    Top = 67
    Width = 8
    Height = 13
    Caption = 'm'
  end
  object Label6: TLabel
    Left = 16
    Top = 99
    Width = 6
    Height = 13
    Caption = 'Z'
  end
  object Label7: TLabel
    Left = 160
    Top = 94
    Width = 8
    Height = 13
    Caption = 'm'
  end
  object EditX: TEdit
    Left = 28
    Top = 37
    Width = 121
    Height = 21
    TabOrder = 0
  end
  object EditY: TEdit
    Left = 28
    Top = 64
    Width = 121
    Height = 21
    TabOrder = 1
  end
  object EditZ: TEdit
    Left = 28
    Top = 91
    Width = 121
    Height = 21
    TabOrder = 2
  end
  object Button1: TButton
    Left = 117
    Top = 118
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 3
    OnClick = Button1Click
  end
end
