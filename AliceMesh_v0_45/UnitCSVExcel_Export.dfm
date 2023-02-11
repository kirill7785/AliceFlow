object FormCSV_Excel_export: TFormCSV_Excel_export
  Left = 0
  Top = 0
  Caption = 'CSV/Excel export'
  ClientHeight = 181
  ClientWidth = 258
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object RadioGroup1: TRadioGroup
    Left = 8
    Top = 16
    Width = 233
    Height = 105
    Caption = 'Geometry type for CVS\Excel export'
    ItemIndex = 0
    Items.Strings = (
      'Prism'
      'Cylinder')
    TabOrder = 0
  end
  object ButtonExport: TButton
    Left = 136
    Top = 144
    Width = 75
    Height = 25
    Caption = 'Export'
    TabOrder = 1
    OnClick = ButtonExportClick
  end
end
