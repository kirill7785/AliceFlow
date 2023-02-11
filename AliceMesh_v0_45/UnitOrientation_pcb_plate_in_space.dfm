object FormOrientation_pcb_plate: TFormOrientation_pcb_plate
  Left = 0
  Top = 0
  Caption = 'Orientation pcb plate in space'
  ClientHeight = 132
  ClientWidth = 297
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object RadioGroupOrientation_pcb: TRadioGroup
    Left = 16
    Top = 16
    Width = 185
    Height = 105
    Caption = 'Orientation pcb plate'
    ItemIndex = 0
    Items.Strings = (
      'XOY'
      'XOY->XOZ')
    TabOrder = 0
  end
end
