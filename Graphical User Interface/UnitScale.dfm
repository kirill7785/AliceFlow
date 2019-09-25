object FormScale: TFormScale
  Left = 0
  Top = 0
  Caption = 'Scale'
  ClientHeight = 140
  ClientWidth = 229
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 0
    Top = 0
    Width = 225
    Height = 137
    Caption = 'Scale all geometry object'
    Color = clMoneyGreen
    ParentBackground = False
    ParentColor = False
    TabOrder = 0
    object Label1: TLabel
      Left = 16
      Top = 40
      Width = 46
      Height = 13
      Caption = 'maschtab'
    end
    object Edit1: TEdit
      Left = 80
      Top = 32
      Width = 121
      Height = 21
      TabOrder = 0
      Text = '1'
    end
    object ButtonScale: TButton
      Left = 126
      Top = 88
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = ButtonScaleClick
    end
  end
end
