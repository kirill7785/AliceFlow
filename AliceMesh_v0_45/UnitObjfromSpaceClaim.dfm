object FormObjfromSpaceClaim: TFormObjfromSpaceClaim
  Left = 0
  Top = 0
  Caption = 'Read Obj files from SpaceClaim'
  ClientHeight = 234
  ClientWidth = 311
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object ButtonApply: TButton
    Left = 208
    Top = 200
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 0
    OnClick = ButtonApplyClick
  end
  object GroupBoxextrudesource2D23D: TGroupBox
    Left = 8
    Top = 8
    Width = 217
    Height = 129
    Caption = 'Extrude 2D source to 3D source'
    TabOrder = 1
    object LabelX: TLabel
      Left = 16
      Top = 32
      Width = 38
      Height = 13
      Caption = 'normalX'
    end
    object Label1: TLabel
      Left = 16
      Top = 64
      Width = 38
      Height = 13
      Caption = 'normalY'
    end
    object Label2: TLabel
      Left = 16
      Top = 96
      Width = 38
      Height = 13
      Caption = 'normalZ'
    end
    object ComboBoxnormalX: TComboBox
      Left = 68
      Top = 29
      Width = 125
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = 'thickness +25%'
      Items.Strings = (
        'thickness +25%'
        'thickness  -25%')
    end
    object ComboBoxnormalY: TComboBox
      Left = 68
      Top = 56
      Width = 125
      Height = 21
      ItemIndex = 0
      TabOrder = 1
      Text = 'thickness +25%'
      Items.Strings = (
        'thickness +25%'
        'thickness  -25%')
    end
    object ComboBoxnormalZ: TComboBox
      Left = 68
      Top = 93
      Width = 125
      Height = 21
      ItemIndex = 0
      TabOrder = 2
      Text = 'thickness +25%'
      Items.Strings = (
        'thickness +25%'
        'thickness  -25%')
    end
  end
end
