object FormPattern: TFormPattern
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Text name fragment Pattern'
  ClientHeight = 186
  ClientWidth = 393
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
    Left = 8
    Top = 128
    Width = 135
    Height = 13
    Caption = 'text name fragment pattern'
  end
  object RadioGroupPattern: TRadioGroup
    Left = 0
    Top = 0
    Width = 233
    Height = 105
    Caption = 'RadioGroupPattern'
    ItemIndex = 0
    Items.Strings = (
      'Only Select Item'
      'On include text name fragment pattern')
    TabOrder = 0
  end
  object EdittextnamefragmentPattern: TEdit
    Left = 168
    Top = 125
    Width = 225
    Height = 21
    TabOrder = 1
  end
  object ButtonApply: TButton
    Left = 288
    Top = 161
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 2
    OnClick = ButtonApplyClick
  end
end
