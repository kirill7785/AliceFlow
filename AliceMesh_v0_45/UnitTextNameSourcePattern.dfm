object FormTextNameSourcePattern: TFormTextNameSourcePattern
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Set Text pattern for source'
  ClientHeight = 80
  ClientWidth = 369
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
    Top = 0
    Width = 200
    Height = 13
    Caption = 'Enter Text pattern for heat power source'
  end
  object EditSourcepatternname: TEdit
    Left = 0
    Top = 19
    Width = 369
    Height = 21
    TabOrder = 0
  end
  object ButtonApply: TButton
    Left = 264
    Top = 55
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 1
    OnClick = ButtonApplyClick
  end
end
