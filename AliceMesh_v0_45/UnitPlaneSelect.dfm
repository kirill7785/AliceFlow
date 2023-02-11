object FormSelectPlaneRotation: TFormSelectPlaneRotation
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Select Plane rotation'
  ClientHeight = 25
  ClientWidth = 185
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
    Left = 0
    Top = 8
    Width = 26
    Height = 13
    Caption = 'Plane'
  end
  object ComboBoxPlane: TComboBox
    Left = 40
    Top = 2
    Width = 49
    Height = 21
    ItemIndex = 0
    TabOrder = 0
    Text = 'XY'
    Items.Strings = (
      'XY'
      'XZ'
      'YZ')
  end
  object ButtonApply: TButton
    Left = 110
    Top = 0
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 1
    OnClick = ButtonApplyClick
  end
end
