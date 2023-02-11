object Formwallgeometryposition_init: TFormwallgeometryposition_init
  Left = 0
  Top = 0
  Caption = 'Wall quick gemetry position initialization'
  ClientHeight = 185
  ClientWidth = 505
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object RadioGroupwallinitpos: TRadioGroup
    Left = 0
    Top = 0
    Width = 121
    Height = 185
    Caption = 'Wall initial position'
    ItemIndex = 0
    Items.Strings = (
      'None'
      'cabinet min X'
      'cabinet max X'
      'cabinet min Y'
      'cabinet max Y'
      'cabinet min Z'
      'cabinet max Z')
    TabOrder = 0
    OnClick = RadioGroupwallinitposClick
  end
  object ButtonNext: TButton
    Left = 127
    Top = 152
    Width = 75
    Height = 25
    Caption = 'Next'
    TabOrder = 1
    OnClick = ButtonNextClick
  end
end
