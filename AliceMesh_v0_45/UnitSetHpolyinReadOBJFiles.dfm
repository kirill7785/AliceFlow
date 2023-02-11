object FormHpolyinOBJFiles: TFormHpolyinOBJFiles
  Left = 0
  Top = 0
  Caption = 'FormHpolyinOBJFiles'
  ClientHeight = 133
  ClientWidth = 345
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
    Top = 24
    Width = 60
    Height = 13
    Caption = 'block name :'
  end
  object Labelreadblockname: TLabel
    Left = 104
    Top = 24
    Width = 97
    Height = 13
    Caption = 'Labelreadblockname'
  end
  object LabeliPlane: TLabel
    Left = 16
    Top = 48
    Width = 102
    Height = 13
    Caption = 'Polygon iPlaneObj2 : '
  end
  object LabelreadiPlaneObj2: TLabel
    Left = 136
    Top = 48
    Width = 98
    Height = 13
    Caption = 'LabelreadiPlaneObj2'
  end
  object Label2: TLabel
    Left = 16
    Top = 80
    Width = 46
    Height = 13
    Caption = 'Set Hpoly'
  end
  object Label3: TLabel
    Left = 216
    Top = 80
    Width = 8
    Height = 13
    Caption = 'm'
  end
  object ButtonApply: TButton
    Left = 262
    Top = 104
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 0
    OnClick = ButtonApplyClick
  end
  object EditHpoly: TEdit
    Left = 80
    Top = 80
    Width = 121
    Height = 21
    TabOrder = 1
    Text = '0'
  end
end
