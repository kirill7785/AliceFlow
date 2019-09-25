object Formamg1r5Parameters: TFormamg1r5Parameters
  Left = 0
  Top = 0
  Caption = 'amg1r5 Parameters'
  ClientHeight = 122
  ClientWidth = 277
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBoxStabilisation: TGroupBox
    Left = 8
    Top = 8
    Width = 273
    Height = 65
    Caption = 'Stabilisation'
    TabOrder = 0
    object ComboBoxStabilization: TComboBox
      Left = 3
      Top = 24
      Width = 198
      Height = 21
      ItemIndex = 1
      TabOrder = 0
      Text = 'BiCGStab'
      Items.Strings = (
        'none'
        'BiCGStab'
        'FGMRes')
    end
  end
  object ButtonApply: TButton
    Left = 194
    Top = 88
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 1
    OnClick = ButtonApplyClick
  end
  object CheckBox_amg1r6: TCheckBox
    Left = 19
    Top = 92
    Width = 78
    Height = 17
    Caption = 'amg1r6'
    TabOrder = 2
  end
end
